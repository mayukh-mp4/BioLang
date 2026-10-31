#include "interpreter.hpp"
#include <algorithm>
#include <cctype>
#include <chrono>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <sstream>

// ---------- sequence operations ----------
// Standard genetic code, bases ordered T, C, A, G (index = 16*b1 + 4*b2 + b3).
static const char* CODON_TABLE =
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

static int baseIndex(char c) {
    switch (c) {
        case 'T': case 'U': return 0;
        case 'C': return 1;
        case 'A': return 2;
        case 'G': return 3;
    }
    return -1;
}

// The input is treated as the coding (sense) strand, so transcription is T -> U.
static std::string transcribe(const std::string& dna) {
    std::string r = dna;
    std::replace(r.begin(), r.end(), 'T', 'U');
    return r;
}

// Reads codons from position 0. Stop codons are kept as '*' (Biopython-compatible);
// a trailing incomplete codon is ignored.
static std::string translate(const std::string& rna) {
    std::string p;
    for (size_t i = 0; i + 2 < rna.size(); i += 3) {
        int a = baseIndex(rna[i]), b = baseIndex(rna[i + 1]), c = baseIndex(rna[i + 2]);
        if (a < 0 || b < 0 || c < 0) { p += 'X'; continue; }  // ambiguous base (e.g. N)
        p += CODON_TABLE[16 * a + 4 * b + c];
    }
    return p;
}

static std::string revcomp(const std::string& dna) {
    std::string r(dna.rbegin(), dna.rend());
    for (char& c : r) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'R': c = 'Y'; break;
            case 'Y': c = 'R'; break;
            case 'K': c = 'M'; break;
            case 'M': c = 'K'; break;
            case 'B': c = 'V'; break;
            case 'V': c = 'B'; break;
            case 'D': c = 'H'; break;
            case 'H': c = 'D'; break;
        }  // S, W and N are their own complements
    }
    return r;
}

static double gcContent(const std::string& s) {
    if (s.empty()) return 0.0;
    size_t gc = 0;
    for (char c : s) if (c == 'G' || c == 'C') gc++;
    return static_cast<double>(gc) / s.size();
}


// Reads the FIRST record of a FASTA file. Letters are upper-cased; IUPAC ambiguity
// codes (N, R, Y, ...) are kept, anything else is a runtime error.
static std::string loadFasta(const std::string& path, const Node& at) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw RuntimeError{"cannot open FASTA file '" + path + "'", at.line, at.col};
    static const std::string allowed = "ACGTNRYKMSWBDHV";
    std::string line, seq;
    bool inRecord = false;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == ';') continue;
        if (line[0] == '>') {
            if (inRecord) break;  // second record: stop
            inRecord = true;
            continue;
        }
        inRecord = true;
        for (char ch : line) {
            if (ch == ' ') continue;
            char u = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            if (allowed.find(u) == std::string::npos)
                throw RuntimeError{std::string("invalid character '") + ch + "' in FASTA file '" + path + "'",
                                   at.line, at.col};
            seq += u;
        }
    }
    if (seq.empty()) throw RuntimeError{"no sequence found in '" + path + "'", at.line, at.col};
    return seq;
}

// ---------- values ----------
static Value makeSeq(Type t, std::string s) { Value v; v.type = t; v.s = std::move(s); return v; }
static Value makeInt(long long i)  { Value v; v.type = Type::INT;   v.i = i; return v; }
static Value makeFloat(double f)   { Value v; v.type = Type::FLOAT; v.f = f; return v; }
static Value makeBool(bool b)      { Value v; v.type = Type::BOOL;  v.b = b; return v; }
static bool isNum(const Value& v)  { return v.type == Type::INT || v.type == Type::FLOAT; }
static double toD(const Value& v)  { return v.type == Type::INT ? static_cast<double>(v.i) : v.f; }

std::string Interpreter::show(const Value& v) {
    std::ostringstream os;
    switch (v.type) {
        case Type::INT:   os << v.i; break;
        case Type::FLOAT: os << v.f; break;
        case Type::BOOL:  os << (v.b ? "true" : "false"); break;
        case Type::MOTIF: os << "/" << v.s << "/"; break;
        case Type::LIST_INT:
            os << "[";
            for (size_t i = 0; i < v.list.size(); i++) os << (i ? ", " : "") << v.list[i];
            os << "]";
            break;
        default: os << v.s; break;
    }
    return os.str();
}

// ---------- interpreter ----------
bool Interpreter::run(const Program& program) {
    env_.clear();
    try {
        for (const auto& s : program.stmts) exec(*s);
    } catch (const RuntimeError& e) {
        std::fprintf(stderr, "runtime error at %d:%d: %s\n", e.line, e.col, e.msg.c_str());
        return false;
    }
    return true;
}

// Gives a value the declared variable type (string literal -> dna, int -> float).
static Value coerce(Value v, Type target) {
    if (target == Type::FLOAT && v.type == Type::INT) return makeFloat(static_cast<double>(v.i));
    if (v.type == Type::STRING &&
        (target == Type::DNA || target == Type::RNA || target == Type::PROTEIN))
        v.type = target;
    return v;
}

static Type declaredType(TokType t) {
    switch (t) {
        case TokType::KW_DNA:     return Type::DNA;
        case TokType::KW_RNA:     return Type::RNA;
        case TokType::KW_PROTEIN: return Type::PROTEIN;
        case TokType::KW_MOTIF:   return Type::MOTIF;
        case TokType::KW_INT:     return Type::INT;
        case TokType::KW_FLOAT:   return Type::FLOAT;
        case TokType::KW_STRING:  return Type::STRING;
        default:                  return Type::BOOL;
    }
}

void Interpreter::exec(const Stmt& s) {
    if (const auto* d = dynamic_cast<const VarDecl*>(&s)) {
        env_[d->name] = coerce(eval(*d->init), declaredType(d->type));
    } else if (const auto* a = dynamic_cast<const Assign*>(&s)) {
        Type t = env_[a->name].type;
        env_[a->name] = coerce(eval(*a->value), t);
    } else if (const auto* p = dynamic_cast<const PrintStmt*>(&s)) {
        std::cout << show(eval(*p->expr)) << "\n";
    } else if (const auto* e = dynamic_cast<const ExprStmt*>(&s)) {
        eval(*e->expr);
    } else if (const auto* i = dynamic_cast<const IfStmt*>(&s)) {
        Value c = eval(*i->cond);
        for (const auto& st : c.b ? i->thenBody : i->elseBody) exec(*st);
    } else if (const auto* f = dynamic_cast<const ForStmt*>(&s)) {
        Value it = eval(*f->iter);
        for (long long v : it.list) {
            env_[f->var] = makeInt(v);
            for (const auto& st : f->body) exec(*st);
        }
    }
}

Value Interpreter::eval(const Expr& e) {
    if (const auto* s = dynamic_cast<const StringLit*>(&e)) { Value v; v.type = Type::STRING; v.s = s->value; return v; }
    if (const auto* r = dynamic_cast<const RegexLit*>(&e)) {
        Value v = makeSeq(Type::MOTIF, r->pattern);
        v.dfa = r->dfa ? r->dfa : compileMotif(r->pattern);
        return v;
    }
    if (const auto* n = dynamic_cast<const NumberLit*>(&e)) {
        if (n->text.find('.') != std::string::npos) return makeFloat(std::stod(n->text));
        return makeInt(std::stoll(n->text));
    }
    if (const auto* id = dynamic_cast<const Ident*>(&e)) return env_.at(id->name);
    if (const auto* b = dynamic_cast<const Binary*>(&e)) return evalBinary(*b);
    if (const auto* ix = dynamic_cast<const Index*>(&e)) {
        Value base = eval(*ix->base), k = eval(*ix->idx);
        if (k.i < 0 || k.i >= static_cast<long long>(base.list.size()))
            throw RuntimeError{"index " + std::to_string(k.i) + " out of range (list has " +
                               std::to_string(base.list.size()) + " items)", ix->line, ix->col};
        return makeInt(base.list[static_cast<size_t>(k.i)]);
    }
    if (const auto* sl = dynamic_cast<const Slice*>(&e)) {
        Value base = eval(*sl->base);
        long long n = static_cast<long long>(base.s.size());
        long long lo = sl->lo ? eval(*sl->lo).i : 0;
        long long hi = sl->hi ? eval(*sl->hi).i : n;
        if (lo < 0 || hi < 0 || lo > n || hi > n || lo > hi)
            throw RuntimeError{"slice [" + std::to_string(lo) + ":" + std::to_string(hi) +
                               "] out of range (length " + std::to_string(n) + ")", sl->line, sl->col};
        return makeSeq(base.type, base.s.substr(static_cast<size_t>(lo), static_cast<size_t>(hi - lo)));
    }
    if (const auto* c = dynamic_cast<const Call*>(&e))   return evalCall(*c);
    throw RuntimeError{"cannot evaluate expression", e.line, e.col};
}

Value Interpreter::evalBinary(const Binary& b) {
    Value l = eval(*b.lhs), r = eval(*b.rhs);
    const std::string& op = b.op;

    if (op == "==" || op == "!=") {
        bool eq = (isNum(l) && isNum(r)) ? toD(l) == toD(r)
                : (l.type == Type::BOOL) ? l.b == r.b
                                         : l.s == r.s;
        return makeBool(op == "==" ? eq : !eq);
    }
    if (op == "<")  return makeBool(toD(l) < toD(r));
    if (op == ">")  return makeBool(toD(l) > toD(r));
    if (op == "<=") return makeBool(toD(l) <= toD(r));
    if (op == ">=") return makeBool(toD(l) >= toD(r));

    if (op == "+" && !isNum(l)) return makeSeq(l.type, l.s + r.s);  // sequence / string concatenation

    if (l.type == Type::INT && r.type == Type::INT) {
        if (op == "+") return makeInt(l.i + r.i);
        if (op == "-") return makeInt(l.i - r.i);
        if (op == "*") return makeInt(l.i * r.i);
        if (r.i == 0) throw RuntimeError{"division by zero", b.line, b.col};
        return makeInt(l.i / r.i);
    }
    double x = toD(l), y = toD(r);
    if (op == "+") return makeFloat(x + y);
    if (op == "-") return makeFloat(x - y);
    if (op == "*") return makeFloat(x * y);
    if (y == 0.0) throw RuntimeError{"division by zero", b.line, b.col};
    return makeFloat(x / y);
}

Value Interpreter::evalCall(const Call& c) {
    const std::string& f = c.callee;

    if (f == "find_all") {
        const auto* in = static_cast<const InExpr*>(c.args[0].get());
        Value m = eval(*in->needle), s = eval(*in->haystack);
        Value out;
        out.type = Type::LIST_INT;
        auto t0 = std::chrono::steady_clock::now();
        out.list = m.dfa->findAll(s.s);
        scanSeconds_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        return out;
    }

    if (f == "match_end") {
        const auto* in = static_cast<const InExpr*>(c.args[0].get());
        Value m = eval(*in->needle), s = eval(*in->haystack), p = eval(*c.args[1]);
        if (p.i < 0 || p.i > static_cast<long long>(s.s.size()))
            throw RuntimeError{"match_end position " + std::to_string(p.i) + " out of range", c.line, c.col};
        return makeInt(m.dfa->matchEnd(s.s, static_cast<size_t>(p.i)));
    }

    Value a = eval(*c.args[0]);
    if (f == "load_fasta") return makeSeq(Type::DNA, loadFasta(a.s, c));
    if (f == "count")      return makeInt(static_cast<long long>(a.list.size()));
    if (f == "transcribe") return makeSeq(Type::RNA, transcribe(a.s));
    if (f == "translate")  return makeSeq(Type::PROTEIN, translate(a.s));
    if (f == "revcomp")    return makeSeq(Type::DNA, revcomp(a.s));
    if (f == "gc_content") return makeFloat(gcContent(a.s));
    if (f == "len")        return makeInt(static_cast<long long>(a.s.size()));
    throw RuntimeError{"unknown function '" + f + "'", c.line, c.col};
}
