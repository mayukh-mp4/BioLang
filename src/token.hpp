#pragma once
#include <string>

enum class TokType {
    // type keywords
    KW_DNA, KW_RNA, KW_PROTEIN, KW_MOTIF, KW_INT, KW_FLOAT, KW_BOOL,
    // other keywords
    KW_IN, KW_IF, KW_ELSE, KW_FOR, KW_PRINT,
    // literals / names
    IDENT, STRING, NUMBER, REGEX,
    // operators and punctuation
    ASSIGN, PLUS, MINUS, STAR, SLASH, EQ, LT, GT,
    LPAREN, RPAREN, COMMA, NEWLINE,
    // special
    END, ERROR
};

struct Token {
    TokType type;
    std::string lexeme;
    int line;
    int col;
};

inline const char* tokName(TokType t) {
    switch (t) {
        case TokType::KW_DNA:     return "KW_DNA";
        case TokType::KW_RNA:     return "KW_RNA";
        case TokType::KW_PROTEIN: return "KW_PROTEIN";
        case TokType::KW_MOTIF:   return "KW_MOTIF";
        case TokType::KW_INT:     return "KW_INT";
        case TokType::KW_FLOAT:   return "KW_FLOAT";
        case TokType::KW_BOOL:    return "KW_BOOL";
        case TokType::KW_IN:      return "KW_IN";
        case TokType::KW_IF:      return "KW_IF";
        case TokType::KW_ELSE:    return "KW_ELSE";
        case TokType::KW_FOR:     return "KW_FOR";
        case TokType::KW_PRINT:   return "KW_PRINT";
        case TokType::IDENT:      return "IDENT";
        case TokType::STRING:     return "STRING";
        case TokType::NUMBER:     return "NUMBER";
        case TokType::REGEX:      return "REGEX";
        case TokType::ASSIGN:     return "ASSIGN";
        case TokType::PLUS:       return "PLUS";
        case TokType::MINUS:      return "MINUS";
        case TokType::STAR:       return "STAR";
        case TokType::SLASH:      return "SLASH";
        case TokType::EQ:         return "EQ";
        case TokType::LT:         return "LT";
        case TokType::GT:         return "GT";
        case TokType::LPAREN:     return "LPAREN";
        case TokType::RPAREN:     return "RPAREN";
        case TokType::COMMA:      return "COMMA";
        case TokType::NEWLINE:    return "NEWLINE";
        case TokType::END:        return "END";
        case TokType::ERROR:      return "ERROR";
    }
    return "?";
}
