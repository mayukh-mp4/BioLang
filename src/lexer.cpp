#include "lexer.hpp"
#include <cctype>
#include <unordered_map>

static const std::unordered_map<std::string, TokType> KEYWORDS = {
    {"dna", TokType::KW_DNA},     {"rna", TokType::KW_RNA},
    {"protein", TokType::KW_PROTEIN}, {"motif", TokType::KW_MOTIF},
    {"int", TokType::KW_INT},     {"float", TokType::KW_FLOAT},
    {"bool", TokType::KW_BOOL},   {"in", TokType::KW_IN},
    {"if", TokType::KW_IF},       {"else", TokType::KW_ELSE},
    {"for", TokType::KW_FOR},     {"print", TokType::KW_PRINT},
};

Lexer::Lexer(std::string source) : src_(std::move(source)) {}

bool Lexer::atEnd() const { return pos_ >= src_.size(); }

char Lexer::peek(size_t offset) const {
    return (pos_ + offset < src_.size()) ? src_[pos_ + offset] : '\0';
}

char Lexer::advance() {
    char c = src_[pos_++];
    if (c == '\n') { line_++; col_ = 1; }
    else           { col_++; }
    return c;
}

Token Lexer::make(TokType t, std::string lexeme, int line, int col) {
    prev_ = t;
    return Token{t, std::move(lexeme), line, col};
}

// Skips spaces, tabs, carriage returns and '#' comments.
// Newlines are NOT skipped: they become NEWLINE tokens.
void Lexer::skipSpaceAndComments() {
    while (!atEnd()) {
        char c = peek();
        if (c == ' ' || c == '\t' || c == '\r') {
            advance();
        } else if (c == '#') {
            while (!atEnd() && peek() != '\n') advance();
        } else {
            break;
        }
    }
}

// A '/' starts a regex literal only where an expression may begin.
// Otherwise it is the division operator.
bool Lexer::regexAllowed() const {
    return prev_ == TokType::ASSIGN || prev_ == TokType::LPAREN ||
           prev_ == TokType::COMMA  || prev_ == TokType::NEWLINE;
}

Token Lexer::readIdent(int line, int col) {
    std::string s;
    while (std::isalnum(static_cast<unsigned char>(peek())) || peek() == '_')
        s += advance();
    auto it = KEYWORDS.find(s);
    return make(it != KEYWORDS.end() ? it->second : TokType::IDENT, s, line, col);
}

Token Lexer::readNumber(int line, int col) {
    std::string s;
    while (std::isdigit(static_cast<unsigned char>(peek()))) s += advance();
    if (peek() == '.' && std::isdigit(static_cast<unsigned char>(peek(1)))) {
        s += advance();  // '.'
        while (std::isdigit(static_cast<unsigned char>(peek()))) s += advance();
    }
    return make(TokType::NUMBER, s, line, col);
}

Token Lexer::readString(int line, int col) {
    advance();  // opening quote
    std::string s;
    while (!atEnd() && peek() != '"' && peek() != '\n') s += advance();
    if (peek() != '"')
        return make(TokType::ERROR, "unterminated string", line, col);
    advance();  // closing quote
    return make(TokType::STRING, s, line, col);
}

Token Lexer::readRegex(int line, int col) {
    advance();  // opening '/'
    std::string s;
    while (!atEnd() && peek() != '/' && peek() != '\n') {
        if (peek() == '\\' && peek(1) != '\0') s += advance();  // keep escape
        s += advance();
    }
    if (peek() != '/')
        return make(TokType::ERROR, "unterminated motif pattern", line, col);
    advance();  // closing '/'
    return make(TokType::REGEX, s, line, col);
}

Token Lexer::next() {
    int line = line_, col = col_;
    char c = peek();

    if (std::isalpha(static_cast<unsigned char>(c)) || c == '_')
        return readIdent(line, col);
    if (std::isdigit(static_cast<unsigned char>(c)))
        return readNumber(line, col);
    if (c == '"')
        return readString(line, col);

    switch (c) {
        case '/':
            if (regexAllowed()) return readRegex(line, col);
            advance();
            return make(TokType::SLASH, "/", line, col);
        case '=':
            advance();
            if (peek() == '=') { advance(); return make(TokType::EQ, "==", line, col); }
            return make(TokType::ASSIGN, "=", line, col);
        case '+': advance(); return make(TokType::PLUS,   "+", line, col);
        case '-': advance(); return make(TokType::MINUS,  "-", line, col);
        case '*': advance(); return make(TokType::STAR,   "*", line, col);
        case '<':
            advance();
            if (peek() == '=') { advance(); return make(TokType::LE, "<=", line, col); }
            return make(TokType::LT, "<", line, col);
        case '>':
            advance();
            if (peek() == '=') { advance(); return make(TokType::GE, ">=", line, col); }
            return make(TokType::GT, ">", line, col);
        case '!':
            if (peek(1) == '=') { advance(); advance(); return make(TokType::NE, "!=", line, col); }
            break;
        case '{': advance(); return make(TokType::LBRACE, "{", line, col);
        case '}': advance(); return make(TokType::RBRACE, "}", line, col);
        case '(': advance(); return make(TokType::LPAREN, "(", line, col);
        case ')': advance(); return make(TokType::RPAREN, ")", line, col);
        case ',': advance(); return make(TokType::COMMA,  ",", line, col);
    }

    advance();
    return make(TokType::ERROR, std::string("unexpected character '") + c + "'", line, col);
}

std::vector<Token> Lexer::tokenize() {
    std::vector<Token> out;
    while (true) {
        skipSpaceAndComments();
        if (atEnd()) break;

        if (peek() == '\n') {
            int line = line_, col = col_;
            advance();
            if (prev_ != TokType::NEWLINE)  // collapse blank lines
                out.push_back(make(TokType::NEWLINE, "\\n", line, col));
            continue;
        }
        out.push_back(next());
    }
    out.push_back(Token{TokType::END, "", line_, col_});
    return out;
}
