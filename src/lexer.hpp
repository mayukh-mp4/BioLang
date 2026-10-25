#pragma once
#include <string>
#include <vector>
#include "token.hpp"

class Lexer {
public:
    explicit Lexer(std::string source);
    std::vector<Token> tokenize();

private:
    std::string src_;
    size_t pos_ = 0;
    int line_ = 1;
    int col_ = 1;
    TokType prev_ = TokType::NEWLINE;  // type of the last token emitted

    bool atEnd() const;
    char peek(size_t offset = 0) const;
    char advance();
    Token make(TokType t, std::string lexeme, int line, int col);

    void skipSpaceAndComments();
    bool regexAllowed() const;

    Token next();
    Token readIdent(int line, int col);
    Token readNumber(int line, int col);
    Token readString(int line, int col);
    Token readRegex(int line, int col);
};
