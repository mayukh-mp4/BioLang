#pragma once
#include <string>
#include <vector>
#include "ast.hpp"
#include "token.hpp"

struct ParseError {
    std::string msg;
    int line, col;
};

class Parser {
public:
    explicit Parser(std::vector<Token> tokens);
    Program parseProgram();
    const std::vector<std::string>& errors() const { return errors_; }

private:
    std::vector<Token> toks_;
    size_t pos_ = 0;
    std::vector<std::string> errors_;

    const Token& peek(size_t offset = 0) const;
    const Token& advance();
    bool check(TokType t) const;
    bool match(TokType t);
    const Token& expect(TokType t, const char* what);
    [[noreturn]] void fail(const std::string& msg) const;
    void skipNewlines();
    void synchronize();

    StmtPtr parseStatement();
    Block parseBlock();
    ExprPtr parseExpr();
    ExprPtr parseAdditive();
    ExprPtr parseTerm();
    ExprPtr parsePrimary();
    ExprPtr parseCall(const Token& name);
};
