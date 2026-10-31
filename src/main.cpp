#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include "interpreter.hpp"
#include "lexer.hpp"
#include "parser.hpp"
#include "typechecker.hpp"

using Clock = std::chrono::steady_clock;
static double msBetween(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

int main(int argc, char** argv) {
    bool tokensOnly = false, showAst = false, checkOnly = false, showDfa = false, timing = false;
    const char* path = nullptr;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--tokens") == 0) tokensOnly = true;
        else if (std::strcmp(argv[i], "--ast") == 0) showAst = true;
        else if (std::strcmp(argv[i], "--check") == 0) checkOnly = true;
        else if (std::strcmp(argv[i], "--dfa") == 0) showDfa = true;
        else if (std::strcmp(argv[i], "--time") == 0) timing = true;
        else if (std::strcmp(argv[i], "--scan=restart") == 0) scanMode = ScanMode::Restart;
        else if (std::strcmp(argv[i], "--scan=linear") == 0) scanMode = ScanMode::Linear;
        else if (std::strcmp(argv[i], "--scan=auto") == 0) scanMode = ScanMode::Auto;
        else path = argv[i];
    }
    if (!path) {
        std::cerr << "usage: biolang [--tokens | --ast | --check | --dfa | --time | --scan=auto|linear|restart] <file.bl>\n";
        return 1;
    }

    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::cerr << "error: cannot open '" << path << "'\n";
        return 1;
    }
    std::stringstream buf;
    buf << in.rdbuf();

    auto t0 = Clock::now();

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
    auto t1 = Clock::now();

    // 3. Type checking (this also compiles every motif to a DFA)
    TypeChecker checker;
    std::vector<std::string> typeErrors = checker.check(program);
    if (!typeErrors.empty()) {
        for (const auto& e : typeErrors) std::fprintf(stderr, "%s\n", e.c_str());
        std::fprintf(stderr, "\n%zu type error(s)\n", typeErrors.size());
        return 1;
    }
    auto t2 = Clock::now();

    if (showDfa) {
        for (const auto& m : checker.motifs()) m->dump(std::cout);
        if (checker.motifs().empty()) std::cout << "(no motifs in this program)\n";
        return 0;
    }
    if (checkOnly) {
        std::printf("Type check passed (%zu statements)\n", program.stmts.size());
        return 0;
    }

    // 4. Execution
    Interpreter interpreter;
    bool ok = interpreter.run(program);
    auto t3 = Clock::now();

    if (timing)
        std::fprintf(stderr, "[time] parse_ms=%.3f compile_ms=%.3f scan_ms=%.3f total_ms=%.3f\n",
                     msBetween(t0, t1), msBetween(t1, t2), interpreter.scanSeconds() * 1000.0,
                     msBetween(t0, t3));
    return ok ? 0 : 1;
}
