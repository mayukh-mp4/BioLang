#include "parser.hpp"
#include <algorithm>

static bool isTypeKeyword(TokType t) {
    return t == TokType::KW_DNA || t == TokType::KW_RNA || t == TokType::KW_PROTEIN ||
           t == TokType::KW_MOTIF || t == TokType::KW_INT || t == TokType::KW_FLOAT ||
           t == TokType::KW_BOOL;
}

static std::string describe(const Token& t) {
    if (t.type == TokType::END) return "end of file";
    if (t.type == TokType::NEWLINE) return "end of line";
    return "'" + t.lexeme + "'";
}

Parser::Parser(std::vector<Token> tokens) : toks_(std::move(tokens)) {}

const Token& Parser::peek(size_t offset) const {
    return toks_[std::min(pos_ + offset, toks_.size() - 1)];
}

const Token& Parser::advance() {
    const Token& t = toks_[pos_];
    if (pos_ < toks_.size() - 1) pos_++;
    return t;
}

bool Parser::check(TokType t) const { return peek().type == t; }

bool Parser::match(TokType t) {
    if (check(t)) { advance(); return true; }
    return false;
}

void Parser::fail(const std::string& msg) const {
    throw ParseError{msg, peek().line, peek().col};
}

const Token& Parser::expect(TokType t, const char* what) {
    if (!check(t))
        fail(std::string("expected ") + what + ", got " + describe(peek()));
    return advance();
}

void Parser::skipNewlines() { while (match(TokType::NEWLINE)) {} }

// After an error, skip to the end of the current line so parsing can continue.
void Parser::synchronize() {
    while (!check(TokType::END) && !check(TokType::NEWLINE)) advance();
}

Program Parser::parseProgram() {
    Program prog;
    skipNewlines();
    while (!check(TokType::END)) {
        try {
            prog.stmts.push_back(parseStatement());
            if (!check(TokType::NEWLINE) && !check(TokType::END))
                fail("expected end of line, got " + describe(peek()));
        } catch (const ParseError& e) {
            errors_.push_back("error at " + std::to_string(e.line) + ":" +
                              std::to_string(e.col) + ": " + e.msg);
            synchronize();
        }
        skipNewlines();
    }
    return prog;
}

StmtPtr Parser::parseStatement() {
    const Token& first = peek();

    if (isTypeKeyword(first.type)) {
        Token typeTok = advance();
        Token name = expect(TokType::IDENT, "a variable name");
        expect(TokType::ASSIGN, "'='");
        ExprPtr init = parseExpr();
        return std::make_unique<VarDecl>(typeTok.type, typeTok.lexeme, name.lexeme,
                                         std::move(init), typeTok.line, typeTok.col);
    }

    if (first.type == TokType::KW_PRINT) {
        Token kw = advance();
        expect(TokType::LPAREN, "'('");
        ExprPtr e = parseExpr();
        expect(TokType::RPAREN, "')'");
        return std::make_unique<PrintStmt>(std::move(e), kw.line, kw.col);
    }

    if (first.type == TokType::IDENT && peek(1).type == TokType::ASSIGN) {
        Token name = advance();
        advance();  // '='
        ExprPtr v = parseExpr();
        return std::make_unique<Assign>(name.lexeme, std::move(v), name.line, name.col);
    }

    int line = first.line, col = first.col;
    return std::make_unique<ExprStmt>(parseExpr(), line, col);
}

// comparison = additive [ (== | < | >) additive ]
ExprPtr Parser::parseExpr() {
    ExprPtr left = parseAdditive();
    if (check(TokType::EQ) || check(TokType::LT) || check(TokType::GT)) {
        Token op = advance();
        ExprPtr right = parseAdditive();
        return std::make_unique<Binary>(op.lexeme, std::move(left), std::move(right),
                                        op.line, op.col);
    }
    return left;
}

ExprPtr Parser::parseAdditive() {
    ExprPtr left = parseTerm();
    while (check(TokType::PLUS) || check(TokType::MINUS)) {
        Token op = advance();
        ExprPtr right = parseTerm();
        left = std::make_unique<Binary>(op.lexeme, std::move(left), std::move(right),
                                        op.line, op.col);
    }
    return left;
}

ExprPtr Parser::parseTerm() {
    ExprPtr left = parsePrimary();
    while (check(TokType::STAR) || check(TokType::SLASH)) {
        Token op = advance();
        ExprPtr right = parsePrimary();
        left = std::make_unique<Binary>(op.lexeme, std::move(left), std::move(right),
                                        op.line, op.col);
    }
    return left;
}

ExprPtr Parser::parsePrimary() {
    const Token& t = peek();
    switch (t.type) {
        case TokType::STRING: { Token k = advance(); return std::make_unique<StringLit>(k.lexeme, k.line, k.col); }
        case TokType::NUMBER: { Token k = advance(); return std::make_unique<NumberLit>(k.lexeme, k.line, k.col); }
        case TokType::REGEX:  { Token k = advance(); return std::make_unique<RegexLit>(k.lexeme, k.line, k.col); }
        case TokType::IDENT: {
            Token k = advance();
            if (check(TokType::LPAREN)) return parseCall(k);
            return std::make_unique<Ident>(k.lexeme, k.line, k.col);
        }
        case TokType::LPAREN: {
            advance();
            ExprPtr e = parseExpr();
            expect(TokType::RPAREN, "')'");
            return e;
        }
        default:
            fail("expected an expression, got " + describe(t));
    }
}

// call = IDENT "(" [ arg { "," arg } ] ")"     arg = expr [ "in" expr ]
ExprPtr Parser::parseCall(const Token& name) {
    expect(TokType::LPAREN, "'('");
    std::vector<ExprPtr> args;
    if (!check(TokType::RPAREN)) {
        do {
            ExprPtr a = parseExpr();
            if (match(TokType::KW_IN)) {
                ExprPtr hay = parseExpr();
                int l = a->line, c = a->col;
                a = std::make_unique<InExpr>(std::move(a), std::move(hay), l, c);
            }
            args.push_back(std::move(a));
        } while (match(TokType::COMMA));
    }
    expect(TokType::RPAREN, "')'");
    return std::make_unique<Call>(name.lexeme, std::move(args), name.line, name.col);
}
