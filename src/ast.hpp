#pragma once
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "token.hpp"

inline void indent(int depth) { std::cout << std::string(depth * 2, ' '); }

struct Node {
    int line, col;
    Node(int l, int c) : line(l), col(c) {}
    virtual ~Node() = default;
    virtual void print(int depth) const = 0;
};
struct Expr : Node { using Node::Node; };
struct Stmt : Node { using Node::Node; };
using ExprPtr = std::unique_ptr<Expr>;
using StmtPtr = std::unique_ptr<Stmt>;

// ---------- expressions ----------
struct Ident : Expr {
    std::string name;
    Ident(std::string n, int l, int c) : Expr(l, c), name(std::move(n)) {}
    void print(int d) const override { indent(d); std::cout << "Ident " << name << "\n"; }
};

struct StringLit : Expr {
    std::string value;
    StringLit(std::string v, int l, int c) : Expr(l, c), value(std::move(v)) {}
    void print(int d) const override { indent(d); std::cout << "String \"" << value << "\"\n"; }
};

struct NumberLit : Expr {
    std::string text;
    NumberLit(std::string t, int l, int c) : Expr(l, c), text(std::move(t)) {}
    void print(int d) const override { indent(d); std::cout << "Number " << text << "\n"; }
};

struct RegexLit : Expr {
    std::string pattern;
    RegexLit(std::string p, int l, int c) : Expr(l, c), pattern(std::move(p)) {}
    void print(int d) const override { indent(d); std::cout << "Motif /" << pattern << "/\n"; }
};

struct Binary : Expr {
    std::string op;
    ExprPtr lhs, rhs;
    Binary(std::string o, ExprPtr a, ExprPtr b, int l, int c)
        : Expr(l, c), op(std::move(o)), lhs(std::move(a)), rhs(std::move(b)) {}
    void print(int d) const override {
        indent(d); std::cout << "Binary " << op << "\n";
        lhs->print(d + 1); rhs->print(d + 1);
    }
};

// `orf in genome`, only valid as a call argument
struct InExpr : Expr {
    ExprPtr needle, haystack;
    InExpr(ExprPtr n, ExprPtr h, int l, int c)
        : Expr(l, c), needle(std::move(n)), haystack(std::move(h)) {}
    void print(int d) const override {
        indent(d); std::cout << "In\n";
        needle->print(d + 1); haystack->print(d + 1);
    }
};

struct Call : Expr {
    std::string callee;
    std::vector<ExprPtr> args;
    Call(std::string f, std::vector<ExprPtr> a, int l, int c)
        : Expr(l, c), callee(std::move(f)), args(std::move(a)) {}
    void print(int d) const override {
        indent(d); std::cout << "Call " << callee << "\n";
        for (const auto& a : args) a->print(d + 1);
    }
};

// ---------- statements ----------
struct VarDecl : Stmt {
    TokType type;
    std::string typeName, name;
    ExprPtr init;
    VarDecl(TokType t, std::string tn, std::string n, ExprPtr i, int l, int c)
        : Stmt(l, c), type(t), typeName(std::move(tn)), name(std::move(n)), init(std::move(i)) {}
    void print(int d) const override {
        indent(d); std::cout << "VarDecl " << typeName << " " << name << "\n";
        init->print(d + 1);
    }
};

struct Assign : Stmt {
    std::string name;
    ExprPtr value;
    Assign(std::string n, ExprPtr v, int l, int c)
        : Stmt(l, c), name(std::move(n)), value(std::move(v)) {}
    void print(int d) const override {
        indent(d); std::cout << "Assign " << name << "\n";
        value->print(d + 1);
    }
};

struct PrintStmt : Stmt {
    ExprPtr expr;
    PrintStmt(ExprPtr e, int l, int c) : Stmt(l, c), expr(std::move(e)) {}
    void print(int d) const override {
        indent(d); std::cout << "Print\n";
        expr->print(d + 1);
    }
};

struct ExprStmt : Stmt {
    ExprPtr expr;
    ExprStmt(ExprPtr e, int l, int c) : Stmt(l, c), expr(std::move(e)) {}
    void print(int d) const override {
        indent(d); std::cout << "ExprStmt\n";
        expr->print(d + 1);
    }
};

struct Program {
    std::vector<StmtPtr> stmts;
    void print() const {
        std::cout << "Program\n";
        for (const auto& s : stmts) s->print(1);
    }
};
