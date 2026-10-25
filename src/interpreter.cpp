#include "interpreter.hpp"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <regex>
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
        }
    }
    return r;
}

static double gcContent(const std::string& s) {
    if (s.empty()) return 0.0;
    size_t gc = 0;
    for (char c : s) if (c == 'G' || c == 'C') gc++;
    return static_cast<double>(gc) / s.size();
}

// TEMPORARY motif engine: std::regex. This will be replaced by the
// Thompson NFA -> subset construction DFA -> single-pass scanner, and
// kept as the correctness baseline to compare the DFA against.
static std::vector<long long> findAll(const std::string& pattern, const std::string& seq,
                                      const Node& at) {
    std::vector<long long> out;
    try {
        std::regex re(pattern);
        for (std::sregex_iterator it(seq.begin(), seq.end(), re), end; it != end; ++it)
            out.push_back(static_cast<long long>(it->position()));
    } catch (const std::regex_error&) {
        throw RuntimeError{"invalid motif pattern /" + pattern + "/", at.line, at.col};
    }
    return out;
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
    }
}

Value Interpreter::eval(const Expr& e) {
    if (const auto* s = dynamic_cast<const StringLit*>(&e)) { Value v; v.type = Type::STRING; v.s = s->value; return v; }
    if (const auto* r = dynamic_cast<const RegexLit*>(&e))  return makeSeq(Type::MOTIF, r->pattern);
    if (const auto* n = dynamic_cast<const NumberLit*>(&e)) {
        if (n->text.find('.') != std::string::npos) return makeFloat(std::stod(n->text));
        return makeInt(std::stoll(n->text));
    }
    if (const auto* id = dynamic_cast<const Ident*>(&e)) return env_.at(id->name);
    if (const auto* b = dynamic_cast<const Binary*>(&e)) return evalBinary(*b);
    if (const auto* c = dynamic_cast<const Call*>(&e))   return evalCall(*c);
    throw RuntimeError{"cannot evaluate expression", e.line, e.col};
}

Value Interpreter::evalBinary(const Binary& b) {
    Value l = eval(*b.lhs), r = eval(*b.rhs);
    const std::string& op = b.op;

    if (op == "==") {
        if (isNum(l) && isNum(r)) return makeBool(toD(l) == toD(r));
        if (l.type == Type::BOOL) return makeBool(l.b == r.b);
        return makeBool(l.s == r.s);
    }
    if (op == "<") return makeBool(toD(l) < toD(r));
    if (op == ">") return makeBool(toD(l) > toD(r));

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
        out.list = findAll(m.s, s.s, c);
        return out;
    }

    Value a = eval(*c.args[0]);
    if (f == "transcribe") return makeSeq(Type::RNA, transcribe(a.s));
    if (f == "translate")  return makeSeq(Type::PROTEIN, translate(a.s));
    if (f == "revcomp")    return makeSeq(Type::DNA, revcomp(a.s));
    if (f == "gc_content") return makeFloat(gcContent(a.s));
    if (f == "len")        return makeInt(static_cast<long long>(a.s.size()));
    throw RuntimeError{"unknown function '" + f + "'", c.line, c.col};
}
