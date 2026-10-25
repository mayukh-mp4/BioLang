#pragma once
#include <map>
#include <string>
#include <vector>
#include "ast.hpp"

enum class Type {
    DNA, RNA, PROTEIN, MOTIF, INT, FLOAT, BOOL,
    STRING,    // a raw string literal, only convertible to dna/rna/protein
    LIST_INT,  // result of find_all
    UNKNOWN    // result of an erroneous expression; suppresses cascading errors
};

const char* typeName(Type t);

class TypeChecker {
public:
    // Returns the list of errors (empty if the program is well-typed).
    std::vector<std::string> check(const Program& program);

private:
    std::map<std::string, Type> symbols_;
    std::vector<std::string> errors_;

    void error(const Node& at, const std::string& msg);
    void checkStmt(const Stmt& s);
    void checkInit(Type target, const Expr& e, const std::string& varName);
    void checkAlphabet(Type target, const StringLit& lit);
    Type typeOf(const Expr& e);
    Type typeOfBinary(const Binary& b);
    Type typeOfCall(const Call& c);
};
