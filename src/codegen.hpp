#pragma once
#include <string>
#include "ast.hpp"

// Thrown for constructs the LLVM backend does not support.
struct CodegenError {
    std::string msg;
    int line, col;
};

// Emits LLVM IR (textual .ll, needs LLVM >= 15 / opaque pointers) for a type-checked program.
// Motifs must already be compiled (the type checker attaches the DFA to each RegexLit).
std::string generateLLVM(const Program& program);
