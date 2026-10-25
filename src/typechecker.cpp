#include "typechecker.hpp"
#include <set>

// ---------- helpers ----------
const char* typeName(Type t) {
    switch (t) {
        case Type::DNA:      return "dna";
        case Type::RNA:      return "rna";
        case Type::PROTEIN:  return "protein";
        case Type::MOTIF:    return "motif";
        case Type::INT:      return "int";
        case Type::FLOAT:    return "float";
        case Type::BOOL:     return "bool";
        case Type::STRING:   return "string";
        case Type::LIST_INT: return "list<int>";
        case Type::UNKNOWN:  return "?";
    }
    return "?";
}

static bool isSeq(Type t) { return t == Type::DNA || t == Type::RNA || t == Type::PROTEIN; }
static bool isNumeric(Type t) { return t == Type::INT || t == Type::FLOAT; }

static Type typeFromToken(TokType t) {
    switch (t) {
        case TokType::KW_DNA:     return Type::DNA;
        case TokType::KW_RNA:     return Type::RNA;
        case TokType::KW_PROTEIN: return Type::PROTEIN;
        case TokType::KW_MOTIF:   return Type::MOTIF;
        case TokType::KW_INT:     return Type::INT;
        case TokType::KW_FLOAT:   return Type::FLOAT;
        case TokType::KW_BOOL:    return Type::BOOL;
        default:                  return Type::UNKNOWN;
    }
}

static std::string alphabetOf(Type t) {
    if (t == Type::DNA) return "ACGT";
    if (t == Type::RNA) return "ACGU";
    return "ACDEFGHIKLMNPQRSTVWY*";  // 20 amino acids + stop
}

static std::string joinTypes(const std::vector<Type>& ts) {
    std::string s;
    for (size_t i = 0; i < ts.size(); i++) {
        if (i) s += " or ";
        s += typeName(ts[i]);
    }
    return s;
}

// ---------- built-in signatures ----------
struct Sig {
    std::vector<Type> accepts;  // allowed types for the single argument
    Type ret;
    const char* hint;           // shown when a DNA/RNA mix-up is detected
};

static const std::map<std::string, Sig>& builtins() {
    static const std::map<std::string, Sig> table = {
        {"transcribe", {{Type::DNA}, Type::RNA, "transcribe() converts dna to rna"}},
        {"translate",  {{Type::RNA}, Type::PROTEIN, "use transcribe() on the dna first"}},
        {"revcomp",    {{Type::DNA}, Type::DNA, ""}},
        {"gc_content", {{Type::DNA, Type::RNA}, Type::FLOAT, ""}},
        {"len",        {{Type::DNA, Type::RNA, Type::PROTEIN}, Type::INT, ""}},
    };
    return table;
}

// ---------- checker ----------
void TypeChecker::error(const Node& at, const std::string& msg) {
    errors_.push_back("error at " + std::to_string(at.line) + ":" +
                      std::to_string(at.col) + ": " + msg);
}

std::vector<std::string> TypeChecker::check(const Program& program) {
    symbols_.clear();
    errors_.clear();
    for (const auto& s : program.stmts) checkStmt(*s);
    return errors_;
}

void TypeChecker::checkAlphabet(Type target, const StringLit& lit) {
    std::string allowed = alphabetOf(target);
    for (char c : lit.value) {
        if (allowed.find(c) == std::string::npos) {
            std::string list;
            for (char a : allowed) { list += a; list += ' '; }
            list.pop_back();
            error(lit, std::string("invalid character '") + c + "' in " +
                       typeName(target) + " literal (allowed: " + list + ")");
            return;
        }
    }
}

// Checks that `e` can initialise / be assigned to a variable of type `target`.
void TypeChecker::checkInit(Type target, const Expr& e, const std::string& varName) {
    Type t = typeOf(e);
    if (t == Type::UNKNOWN || t == target) return;
    if (target == Type::FLOAT && t == Type::INT) return;  // widening
    if (t == Type::STRING && isSeq(target)) {
        if (const auto* lit = dynamic_cast<const StringLit*>(&e)) checkAlphabet(target, *lit);
        return;
    }
    error(e, std::string("cannot assign ") + typeName(t) + " to " + typeName(target) +
                 " variable '" + varName + "'");
}

void TypeChecker::checkStmt(const Stmt& s) {
    if (const auto* d = dynamic_cast<const VarDecl*>(&s)) {
        Type declared = typeFromToken(d->type);
        checkInit(declared, *d->init, d->name);  // init is checked before the name exists
        if (symbols_.count(d->name))
            error(*d, "variable '" + d->name + "' is already declared");
        else
            symbols_[d->name] = declared;
    } else if (const auto* a = dynamic_cast<const Assign*>(&s)) {
        auto it = symbols_.find(a->name);
        if (it == symbols_.end()) {
            error(*a, "assignment to undeclared variable '" + a->name + "'");
            typeOf(*a->value);
        } else {
            checkInit(it->second, *a->value, a->name);
        }
    } else if (const auto* p = dynamic_cast<const PrintStmt*>(&s)) {
        typeOf(*p->expr);
    } else if (const auto* e = dynamic_cast<const ExprStmt*>(&s)) {
        typeOf(*e->expr);
    }
}

Type TypeChecker::typeOf(const Expr& e) {
    if (dynamic_cast<const StringLit*>(&e)) return Type::STRING;
    if (dynamic_cast<const RegexLit*>(&e))  return Type::MOTIF;
    if (const auto* n = dynamic_cast<const NumberLit*>(&e))
        return n->text.find('.') != std::string::npos ? Type::FLOAT : Type::INT;
    if (const auto* id = dynamic_cast<const Ident*>(&e)) {
        auto it = symbols_.find(id->name);
        if (it == symbols_.end()) {
            error(*id, "undefined variable '" + id->name + "'");
            return Type::UNKNOWN;
        }
        return it->second;
    }
    if (const auto* b = dynamic_cast<const Binary*>(&e)) return typeOfBinary(*b);
    if (const auto* c = dynamic_cast<const Call*>(&e))   return typeOfCall(*c);
    if (const auto* in = dynamic_cast<const InExpr*>(&e)) {
        typeOf(*in->needle);
        typeOf(*in->haystack);
        error(*in, "'in' can only be used inside find_all()");
        return Type::UNKNOWN;
    }
    return Type::UNKNOWN;
}

Type TypeChecker::typeOfBinary(const Binary& b) {
    Type l = typeOf(*b.lhs), r = typeOf(*b.rhs);
    if (l == Type::UNKNOWN || r == Type::UNKNOWN) return Type::UNKNOWN;

    const std::string& op = b.op;
    std::string pair = std::string(typeName(l)) + " and " + typeName(r);

    if (op == "+" || op == "-" || op == "*" || op == "/") {
        if (isNumeric(l) && isNumeric(r))
            return (l == Type::FLOAT || r == Type::FLOAT) ? Type::FLOAT : Type::INT;
        error(b, "operator '" + op + "' cannot be applied to " + pair);
        return Type::UNKNOWN;
    }
    if (op == "==") {
        if (l == r || (isNumeric(l) && isNumeric(r))) return Type::BOOL;
        error(b, std::string("cannot compare ") + typeName(l) + " with " + typeName(r) +
                     ": incompatible types");
        return Type::UNKNOWN;
    }
    // '<' and '>'
    if (isNumeric(l) && isNumeric(r)) return Type::BOOL;
    error(b, "operator '" + op + "' cannot be applied to " + pair);
    return Type::UNKNOWN;
}

Type TypeChecker::typeOfCall(const Call& c) {
    // find_all(motif in sequence) has its own argument form
    if (c.callee == "find_all") {
        const InExpr* in = c.args.size() == 1 ? dynamic_cast<const InExpr*>(c.args[0].get()) : nullptr;
        if (!in) {
            for (const auto& a : c.args) typeOf(*a);
            error(c, "find_all() expects one argument of the form: motif in sequence");
            return Type::LIST_INT;
        }
        Type m = typeOf(*in->needle), s = typeOf(*in->haystack);
        if (m != Type::UNKNOWN && m != Type::MOTIF)
            error(*in->needle, std::string("find_all() expects a motif before 'in', got ") + typeName(m));
        if (s != Type::UNKNOWN && !isSeq(s))
            error(*in->haystack, std::string("find_all() expects dna, rna or protein after 'in', got ") + typeName(s));
        return Type::LIST_INT;
    }

    auto it = builtins().find(c.callee);
    if (it == builtins().end()) {
        for (const auto& a : c.args) typeOf(*a);
        error(c, "unknown function '" + c.callee + "'");
        return Type::UNKNOWN;
    }
    const Sig& sig = it->second;

    std::vector<Type> argTypes;
    for (const auto& a : c.args) argTypes.push_back(typeOf(*a));

    if (c.args.size() != 1) {
        error(c, c.callee + "() takes 1 argument, got " + std::to_string(c.args.size()));
        return sig.ret;
    }
    Type at = argTypes[0];
    if (at == Type::UNKNOWN) return sig.ret;

    for (Type ok : sig.accepts)
        if (at == ok) return sig.ret;

    std::string msg = c.callee + "() expects " + joinTypes(sig.accepts) + ", got " + typeName(at);
    if (at == Type::STRING) msg += " literal (assign it to a typed variable first)";
    else if (*sig.hint && isSeq(at)) msg += "  [hint: " + std::string(sig.hint) + "]";
    error(*c.args[0], msg);
    return sig.ret;  // keep going with the expected type to avoid cascading errors
}
