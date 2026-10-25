#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include "lexer.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: biolang <file.bl>\n";
        return 1;
    }

    std::ifstream in(argv[1], std::ios::binary);
    if (!in) {
        std::cerr << "error: cannot open '" << argv[1] << "'\n";
        return 1;
    }
    std::stringstream buf;
    buf << in.rdbuf();

    Lexer lexer(buf.str());
    std::vector<Token> tokens = lexer.tokenize();

    int errors = 0;
    for (const Token& t : tokens) {
        std::printf("%3d:%-3d  %-10s '%s'\n", t.line, t.col, tokName(t.type),
                    t.lexeme.c_str());
        if (t.type == TokType::ERROR) errors++;
    }

    if (errors) {
        std::fprintf(stderr, "\n%d lexical error(s)\n", errors);
        return 1;
    }
    return 0;
}
