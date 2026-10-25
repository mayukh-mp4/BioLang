#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include "lexer.hpp"
#include "parser.hpp"
#include "typechecker.hpp"

int main(int argc, char** argv) {
    bool tokensOnly = false, showAst = false;
    const char* path = nullptr;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--tokens") == 0) tokensOnly = true;
        else if (std::strcmp(argv[i], "--ast") == 0) showAst = true;
        else path = argv[i];
    }
    if (!path) {
        std::cerr << "usage: biolang [--tokens | --ast] <file.bl>\n";
        return 1;
    }

    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::cerr << "error: cannot open '" << path << "'\n";
        return 1;
    }
    std::stringstream buf;
    buf << in.rdbuf();

    // 1. Lexing
    Lexer lexer(buf.str());
    std::vector<Token> tokens = lexer.tokenize();

    int lexErrors = 0;
    for (const Token& t : tokens) {
        if (tokensOnly)
            std::printf("%3d:%-3d  %-10s '%s'\n", t.line, t.col, tokName(t.type), t.lexeme.c_str());
        if (t.type == TokType::ERROR) {
            std::fprintf(stderr, "lexical error at %d:%d: %s\n", t.line, t.col, t.lexeme.c_str());
            lexErrors++;
        }
    }
    if (lexErrors) return 1;
    if (tokensOnly) return 0;

    // 2. Parsing
    Parser parser(std::move(tokens));
    Program program = parser.parseProgram();
    if (!parser.errors().empty()) {
        for (const auto& e : parser.errors()) std::fprintf(stderr, "%s\n", e.c_str());
        return 1;
    }
    if (showAst) { program.print(); return 0; }

    // 3. Type checking
    TypeChecker checker;
    std::vector<std::string> typeErrors = checker.check(program);
    if (!typeErrors.empty()) {
        for (const auto& e : typeErrors) std::fprintf(stderr, "%s\n", e.c_str());
        std::fprintf(stderr, "\n%zu type error(s)\n", typeErrors.size());
        return 1;
    }
    std::printf("Type check passed (%zu statements)\n", program.stmts.size());
    return 0;
}
