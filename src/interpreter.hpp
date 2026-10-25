#pragma once
#include <map>
#include <string>
#include <vector>
#include "ast.hpp"
#include "motif.hpp"
#include "typechecker.hpp"  // for Type

struct Value {
    Type type = Type::UNKNOWN;
    std::string s;               // dna / rna / protein / string / motif pattern
    long long i = 0;             // int
    double f = 0.0;              // float
    bool b = false;              // bool
    std::vector<long long> list; // list<int>
    std::shared_ptr<Dfa> dfa;    // compiled motif
};

struct RuntimeError {
    std::string msg;
    int line, col;
};

class Interpreter {
public:
    // Runs a type-checked program. Returns false on a runtime error.
    bool run(const Program& program);

    // Total time spent inside find_all scans (for --time / benchmarks).
    double scanSeconds() const { return scanSeconds_; }

private:
    std::map<std::string, Value> env_;
    double scanSeconds_ = 0.0;

    void exec(const Stmt& s);
    Value eval(const Expr& e);
    Value evalBinary(const Binary& b);
    Value evalCall(const Call& c);
    static std::string show(const Value& v);
};
