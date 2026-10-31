// Unit tests for the BioLang lexer, parser and type checker.
// No external framework: a tiny TEST/CHECK harness. Links against all of src/ except main.cpp.
#include <cstdio>
#include <functional>
#include <string>
#include <vector>
#include "lexer.hpp"
#include "parser.hpp"
#include "typechecker.hpp"

// ---------- mini harness ----------
struct TestCase { const char* name; std::function<void()> fn; };
static std::vector<TestCase>& registry() { static std::vector<TestCase> r; return r; }
static int g_failedChecks = 0;
static const char* g_current = "";

struct Registrar {
    Registrar(const char* n, std::function<void()> f) { registry().push_back({n, std::move(f)}); }
};
#define TEST(name) \
    static void test_##name(); \
    static Registrar reg_##name(#name, test_##name); \
    static void test_##name()
#define CHECK(cond) do { if (!(cond)) { \
    std::printf("  FAIL [%s] line %d: %s\n", g_current, __LINE__, #cond); g_failedChecks++; } } while (0)
#define CHECK_EQ(a, b) do { auto _a = (a); auto _b = (b); if (!(_a == _b)) { \
    std::printf("  FAIL [%s] line %d: %s == %s\n", g_current, __LINE__, #a, #b); g_failedChecks++; } } while (0)

// ---------- helpers ----------
static std::vector<Token> lex(const std::string& src) { return Lexer(src).tokenize(); }

static std::vector<TokType> types(const std::string& src) {
    std::vector<TokType> out;
    for (const Token& t : lex(src)) out.push_back(t.type);
    return out;
}

struct Parsed { Program prog; std::vector<std::string> errors; };
static Parsed parse(const std::string& src) {
    Parser p(lex(src));
    Parsed r;
    r.prog = p.parseProgram();
    r.errors = p.errors();
    return r;
}

static std::vector<std::string> typeErrors(const std::string& src) {
    Parsed p = parse(src);
    if (!p.errors.empty()) return {"PARSE ERROR: " + p.errors[0]};
    TypeChecker tc;
    return tc.check(p.prog);
}

static bool anyHas(const std::vector<std::string>& v, const std::string& sub) {
    for (const auto& s : v) if (s.find(sub) != std::string::npos) return true;
    return false;
}

template <class T, class U> static const T* as(const U* p) { return dynamic_cast<const T*>(p); }

// ======================= LEXER =======================
TEST(lexer_keywords_and_identifiers) {
    auto t = lex("dna rna protein motif int float bool string in if else for print foo _bar x1");
    CHECK_EQ(t[0].type, TokType::KW_DNA);
    CHECK_EQ(t[1].type, TokType::KW_RNA);
    CHECK_EQ(t[2].type, TokType::KW_PROTEIN);
    CHECK_EQ(t[3].type, TokType::KW_MOTIF);
    CHECK_EQ(t[4].type, TokType::KW_INT);
    CHECK_EQ(t[5].type, TokType::KW_FLOAT);
    CHECK_EQ(t[6].type, TokType::KW_BOOL);
    CHECK_EQ(t[7].type, TokType::KW_STRING);
    CHECK_EQ(t[8].type, TokType::KW_IN);
    CHECK_EQ(t[9].type, TokType::KW_IF);
    CHECK_EQ(t[10].type, TokType::KW_ELSE);
    CHECK_EQ(t[11].type, TokType::KW_FOR);
    CHECK_EQ(t[12].type, TokType::KW_PRINT);
    CHECK_EQ(t[13].type, TokType::IDENT);
    CHECK_EQ(t[14].lexeme, std::string("_bar"));
    CHECK_EQ(t[15].lexeme, std::string("x1"));
    CHECK_EQ(t.back().type, TokType::END);
}

TEST(lexer_keyword_prefix_is_identifier) {
    auto t = lex("dnax integer iff");
    CHECK_EQ(t[0].type, TokType::IDENT);
    CHECK_EQ(t[1].type, TokType::IDENT);
    CHECK_EQ(t[2].type, TokType::IDENT);
}

TEST(lexer_numbers) {
    auto t = lex("42 3.14 7.");
    CHECK_EQ(t[0].type, TokType::NUMBER);
    CHECK_EQ(t[0].lexeme, std::string("42"));
    CHECK_EQ(t[1].lexeme, std::string("3.14"));
    CHECK_EQ(t[2].lexeme, std::string("7"));          // "7." : the dot needs a digit after it
    CHECK_EQ(t[3].type, TokType::ERROR);              // lone '.' is not a token
}

TEST(lexer_string_literal) {
    auto t = lex("\"ACGT\"");
    CHECK_EQ(t[0].type, TokType::STRING);
    CHECK_EQ(t[0].lexeme, std::string("ACGT"));
}

TEST(lexer_unterminated_string) {
    auto t = lex("dna d = \"ACGT\nprint(d)");
    bool found = false;
    for (const Token& k : t) if (k.type == TokType::ERROR && k.lexeme == "unterminated string") found = true;
    CHECK(found);
}

TEST(lexer_regex_after_assign_paren_comma) {
    CHECK_EQ(lex("m = /AT+/")[2].type, TokType::REGEX);
    CHECK_EQ(lex("find_all(/AT/ in g)")[2].type, TokType::REGEX);
    CHECK_EQ(lex("f(x, /AT/)")[4].type, TokType::REGEX);
    CHECK_EQ(lex("/AT/")[0].type, TokType::REGEX);            // start of file
}

TEST(lexer_slash_is_division_after_operand) {
    auto t = lex("x = a / b");
    CHECK_EQ(t[3].type, TokType::SLASH);
    t = lex("x = 4 / 2");
    CHECK_EQ(t[3].type, TokType::SLASH);
}

// Documents current behaviour: a regex literal is NOT recognised right after '=='.
TEST(lexer_slash_after_equality_is_division_known_quirk) {
    auto t = lex("m == /A/");
    CHECK_EQ(t[2].type, TokType::SLASH);
}

TEST(lexer_regex_keeps_escapes_and_content) {
    auto t = lex("m = /a\\/b[ACG]{2,3}/");
    CHECK_EQ(t[2].type, TokType::REGEX);
    CHECK_EQ(t[2].lexeme, std::string("a\\/b[ACG]{2,3}"));
}

TEST(lexer_regex_with_codon_class_angle_brackets) {
    auto t = lex("m = /ATG(<sense>){10,}<stop>/");
    CHECK_EQ(t[2].type, TokType::REGEX);
    CHECK_EQ(t[2].lexeme, std::string("ATG(<sense>){10,}<stop>"));
}

TEST(lexer_unterminated_regex) {
    auto t = lex("m = /ATG");
    CHECK_EQ(t[2].type, TokType::ERROR);
    CHECK_EQ(t[2].lexeme, std::string("unterminated motif pattern"));
}

TEST(lexer_comparison_operators) {
    auto t = types("a == b != c <= d >= e < f > g");
    CHECK_EQ(t[1], TokType::EQ);
    CHECK_EQ(t[3], TokType::NE);
    CHECK_EQ(t[5], TokType::LE);
    CHECK_EQ(t[7], TokType::GE);
    CHECK_EQ(t[9], TokType::LT);
    CHECK_EQ(t[11], TokType::GT);
}

TEST(lexer_assign_vs_equals) {
    auto t = types("x = y == z");
    CHECK_EQ(t[1], TokType::ASSIGN);
    CHECK_EQ(t[3], TokType::EQ);
}

TEST(lexer_lone_bang_is_error) {
    CHECK_EQ(lex("!x")[0].type, TokType::ERROR);
}

TEST(lexer_unexpected_character) {
    auto t = lex("x = @");
    CHECK_EQ(t[2].type, TokType::ERROR);
    CHECK(t[2].lexeme.find("'@'") != std::string::npos);
}

TEST(lexer_punctuation) {
    auto t = types("{ } ( ) [ ] : ,");
    CHECK_EQ(t[0], TokType::LBRACE);   CHECK_EQ(t[1], TokType::RBRACE);
    CHECK_EQ(t[2], TokType::LPAREN);   CHECK_EQ(t[3], TokType::RPAREN);
    CHECK_EQ(t[4], TokType::LBRACKET); CHECK_EQ(t[5], TokType::RBRACKET);
    CHECK_EQ(t[6], TokType::COLON);    CHECK_EQ(t[7], TokType::COMMA);
}

TEST(lexer_comments_are_skipped) {
    auto t = types("x = 1 # a comment with /slashes/ and \"quotes\"\ny = 2");
    // x = 1 NEWLINE y = 2 END
    CHECK_EQ(t.size(), size_t(8));
    CHECK_EQ(t[3], TokType::NEWLINE);
}

TEST(lexer_blank_lines_collapse_to_one_newline) {
    auto t = types("x = 1\n\n\n\ny = 2");
    int newlines = 0;
    for (TokType k : t) if (k == TokType::NEWLINE) newlines++;
    CHECK_EQ(newlines, 1);
}

TEST(lexer_leading_blank_lines_produce_no_newline_token) {
    auto t = types("\n\nx = 1");
    CHECK_EQ(t[0], TokType::IDENT);
}

TEST(lexer_line_and_column_tracking) {
    auto t = lex("x = 1\n  yy = 22");
    CHECK_EQ(t[0].line, 1); CHECK_EQ(t[0].col, 1);
    CHECK_EQ(t[2].col, 5);
    // tokens: x = 1 NEWLINE yy = 22
    CHECK_EQ(t[4].line, 2); CHECK_EQ(t[4].col, 3);
    CHECK_EQ(t[6].col, 8);
}

TEST(lexer_windows_line_endings) {
    auto t = types("x = 1\r\ny = 2\r\n");   // one NEWLINE between the lines, one trailing
    int newlines = 0;
    for (TokType k : t) if (k == TokType::NEWLINE) newlines++;
    CHECK_EQ(newlines, 2);
    for (TokType k : t) CHECK(k != TokType::ERROR);
}

TEST(lexer_empty_input) {
    auto t = lex("");
    CHECK_EQ(t.size(), size_t(1));
    CHECK_EQ(t[0].type, TokType::END);
}

TEST(lexer_slice_tokens) {
    auto t = types("s[1:3]");
    CHECK_EQ(t[1], TokType::LBRACKET);
    CHECK_EQ(t[3], TokType::COLON);
    CHECK_EQ(t[5], TokType::RBRACKET);
}

// ======================= PARSER =======================
TEST(parser_var_decl) {
    auto p = parse("dna d = \"ACGT\"");
    CHECK(p.errors.empty());
    CHECK_EQ(p.prog.stmts.size(), size_t(1));
    const auto* d = as<VarDecl>(p.prog.stmts[0].get());
    CHECK(d != nullptr);
    if (d) {
        CHECK_EQ(d->name, std::string("d"));
        CHECK_EQ(d->type, TokType::KW_DNA);
        CHECK(as<StringLit>(d->init.get()) != nullptr);
    }
}

TEST(parser_all_declaration_types) {
    auto p = parse("int a = 1\nfloat b = 2.5\nbool c = 1 == 1\nstring s = \"x\"\nmotif m = /A/\nrna r = \"ACGU\"\nprotein q = \"MK\"");
    CHECK(p.errors.empty());
    CHECK_EQ(p.prog.stmts.size(), size_t(7));
}

TEST(parser_assignment_vs_expression_statement) {
    auto p = parse("x = 5\nlen(d)");
    CHECK(p.errors.empty());
    CHECK(as<Assign>(p.prog.stmts[0].get()) != nullptr);
    CHECK(as<ExprStmt>(p.prog.stmts[1].get()) != nullptr);
}

TEST(parser_print) {
    auto p = parse("print(1 + 2)");
    CHECK(p.errors.empty());
    CHECK(as<PrintStmt>(p.prog.stmts[0].get()) != nullptr);
}

TEST(parser_precedence_mul_over_add) {
    auto p = parse("int x = 1 + 2 * 3");
    CHECK(p.errors.empty());
    const auto* d = as<VarDecl>(p.prog.stmts[0].get());
    const auto* top = as<Binary>(d->init.get());
    CHECK(top && top->op == "+");
    if (top) {
        CHECK(as<NumberLit>(top->lhs.get()) != nullptr);
        const auto* rhs = as<Binary>(top->rhs.get());
        CHECK(rhs && rhs->op == "*");
    }
}

TEST(parser_parentheses_override_precedence) {
    auto p = parse("int x = (1 + 2) * 3");
    const auto* d = as<VarDecl>(p.prog.stmts[0].get());
    const auto* top = as<Binary>(d->init.get());
    CHECK(top && top->op == "*");
    if (top) CHECK(as<Binary>(top->lhs.get()) != nullptr);
}

TEST(parser_left_associativity) {
    auto p = parse("int x = 10 - 3 - 2");
    const auto* d = as<VarDecl>(p.prog.stmts[0].get());
    const auto* top = as<Binary>(d->init.get());
    CHECK(top && top->op == "-");
    if (top) {
        CHECK(as<Binary>(top->lhs.get()) != nullptr);   // (10 - 3) - 2
        CHECK(as<NumberLit>(top->rhs.get()) != nullptr);
    }
}

TEST(parser_comparison_binds_looser_than_arithmetic) {
    auto p = parse("bool b = 1 + 2 < 2 * 3");
    CHECK(p.errors.empty());
    const auto* d = as<VarDecl>(p.prog.stmts[0].get());
    const auto* top = as<Binary>(d->init.get());
    CHECK(top && top->op == "<");
}

TEST(parser_chained_comparison_rejected) {
    auto p = parse("bool b = 1 < 2 < 3");
    CHECK(!p.errors.empty());
}

TEST(parser_all_comparison_operators) {
    for (const char* op : {"==", "!=", "<", ">", "<=", ">="}) {
        auto p = parse(std::string("bool b = 1 ") + op + " 2");
        CHECK(p.errors.empty());
    }
}

TEST(parser_slice_forms) {
    for (const char* src : {"s[1:3]", "s[:3]", "s[1:]", "s[:]"}) {
        auto p = parse(src);
        CHECK(p.errors.empty());
        CHECK(!p.prog.stmts.empty() && as<ExprStmt>(p.prog.stmts[0].get()) != nullptr);
        const auto* es = as<ExprStmt>(p.prog.stmts[0].get());
        CHECK(es && as<Slice>(es->expr.get()) != nullptr);
    }
}

TEST(parser_slice_bounds_presence) {
    auto p = parse("s[:3]");
    const auto* sl = as<Slice>(as<ExprStmt>(p.prog.stmts[0].get())->expr.get());
    CHECK(sl && sl->lo == nullptr && sl->hi != nullptr);
    p = parse("s[2:]");
    sl = as<Slice>(as<ExprStmt>(p.prog.stmts[0].get())->expr.get());
    CHECK(sl && sl->lo != nullptr && sl->hi == nullptr);
}

TEST(parser_index_vs_slice) {
    auto p = parse("l[2]");
    CHECK(p.errors.empty());
    CHECK(as<Index>(as<ExprStmt>(p.prog.stmts[0].get())->expr.get()) != nullptr);
}

TEST(parser_empty_brackets_rejected) {
    auto p = parse("s[]");
    CHECK(!p.errors.empty());
}

TEST(parser_slice_with_expression_bounds) {
    auto p = parse("s[i + 1:i + 4]");
    CHECK(p.errors.empty());
}

TEST(parser_call_no_args_and_many_args) {
    auto p = parse("f()\ng(1, 2, 3)");
    CHECK(p.errors.empty());
    const auto* c0 = as<Call>(as<ExprStmt>(p.prog.stmts[0].get())->expr.get());
    const auto* c1 = as<Call>(as<ExprStmt>(p.prog.stmts[1].get())->expr.get());
    CHECK(c0 && c0->args.empty());
    CHECK(c1 && c1->args.size() == 3);
}

TEST(parser_find_all_in_expression) {
    auto p = parse("find_all(/AT/ in g)");
    CHECK(p.errors.empty());
    const auto* c = as<Call>(as<ExprStmt>(p.prog.stmts[0].get())->expr.get());
    CHECK(c && c->callee == "find_all" && c->args.size() == 1);
    if (c && c->args.size() == 1) {
        const auto* in = as<InExpr>(c->args[0].get());
        CHECK(in != nullptr);
        if (in) {
            CHECK(as<RegexLit>(in->needle.get()) != nullptr);
            CHECK(as<Ident>(in->haystack.get()) != nullptr);
        }
    }
}

TEST(parser_match_end_two_args) {
    auto p = parse("match_end(m in g, 3)");
    CHECK(p.errors.empty());
    const auto* c = as<Call>(as<ExprStmt>(p.prog.stmts[0].get())->expr.get());
    CHECK(c && c->args.size() == 2);
    if (c && c->args.size() == 2) {
        CHECK(as<InExpr>(c->args[0].get()) != nullptr);
        CHECK(as<NumberLit>(c->args[1].get()) != nullptr);
    }
}

TEST(parser_if_without_else) {
    auto p = parse("if 1 == 1 {\nprint(1)\n}");
    CHECK(p.errors.empty());
    const auto* i = as<IfStmt>(p.prog.stmts[0].get());
    CHECK(i && i->thenBody.size() == 1 && i->elseBody.empty());
}

TEST(parser_if_else) {
    auto p = parse("if 1 == 1 {\nprint(1)\n} else {\nprint(2)\n}");
    CHECK(p.errors.empty());
    const auto* i = as<IfStmt>(p.prog.stmts[0].get());
    CHECK(i && i->thenBody.size() == 1 && i->elseBody.size() == 1);
}

TEST(parser_else_if_chain_nests) {
    auto p = parse("if 1 == 2 {\nprint(1)\n} else if 1 == 1 {\nprint(2)\n} else {\nprint(3)\n}");
    CHECK(p.errors.empty());
    const auto* i = as<IfStmt>(p.prog.stmts[0].get());
    CHECK(i && i->elseBody.size() == 1);
    if (i && i->elseBody.size() == 1) {
        const auto* inner = as<IfStmt>(i->elseBody[0].get());
        CHECK(inner && inner->elseBody.size() == 1);
    }
}

TEST(parser_for_loop) {
    auto p = parse("for x in find_all(/A/ in g) {\nprint(x)\n}");
    CHECK(p.errors.empty());
    const auto* f = as<ForStmt>(p.prog.stmts[0].get());
    CHECK(f && f->var == "x" && f->body.size() == 1);
}

TEST(parser_closing_brace_on_same_line) {
    auto p = parse("for x in l { print(x) }");
    CHECK(p.errors.empty());
    const auto* f = as<ForStmt>(p.prog.stmts[0].get());
    CHECK(f && f->body.size() == 1);
}

TEST(parser_nested_blocks) {
    auto p = parse("for x in l {\nif x == 1 {\nprint(x)\n}\n}");
    CHECK(p.errors.empty());
    const auto* f = as<ForStmt>(p.prog.stmts[0].get());
    CHECK(f && f->body.size() == 1 && as<IfStmt>(f->body[0].get()) != nullptr);
}

TEST(parser_blank_lines_and_comments_between_statements) {
    auto p = parse("\n# header\n\nint a = 1\n\n# mid\nint b = 2\n\n");
    CHECK(p.errors.empty());
    CHECK_EQ(p.prog.stmts.size(), size_t(2));
}

TEST(parser_error_missing_initializer) {
    auto p = parse("int x =");
    CHECK_EQ(p.errors.size(), size_t(1));
    CHECK(anyHas(p.errors, "expected an expression"));
}

TEST(parser_error_missing_variable_name) {
    auto p = parse("int = 5");
    CHECK(anyHas(p.errors, "expected a variable name"));
}

TEST(parser_error_missing_assign) {
    auto p = parse("int x 5");
    CHECK(anyHas(p.errors, "expected '='"));
}

TEST(parser_error_missing_open_brace) {
    auto p = parse("if 1 == 1\nprint(1)");
    CHECK(anyHas(p.errors, "expected '{'"));
}

TEST(parser_error_unclosed_block) {
    auto p = parse("for x in l {\nprint(x)\n");
    CHECK(anyHas(p.errors, "expected '}'"));
}

TEST(parser_error_unclosed_paren) {
    auto p = parse("print(1 + 2");
    CHECK(anyHas(p.errors, "expected ')'"));
}

TEST(parser_error_trailing_tokens) {
    auto p = parse("int x = 1 2");
    CHECK(anyHas(p.errors, "expected end of line"));
}

TEST(parser_error_position_is_reported) {
    auto p = parse("int a = 1\nint x =");
    CHECK(!p.errors.empty());
    CHECK(p.errors[0].find("error at 2:") == 0);
}

TEST(parser_error_recovery_reports_each_bad_line_and_keeps_good_ones) {
    auto p = parse("int a = \nint b = 2\nint c = \nint d = 4");
    CHECK_EQ(p.errors.size(), size_t(2));
    CHECK_EQ(p.prog.stmts.size(), size_t(2));   // b and d survive
}

TEST(parser_error_loop_variable_must_be_identifier) {
    auto p = parse("for 5 in l {\n}");
    CHECK(anyHas(p.errors, "loop variable"));
}

TEST(parser_error_for_missing_in) {
    auto p = parse("for x l {\n}");
    CHECK(anyHas(p.errors, "expected 'in'"));
}

TEST(parser_always_terminates_on_garbage) {
    // Must not hang or crash.
    auto p = parse("))) ]]] {{{ : , ,");
    CHECK(!p.errors.empty());
}

TEST(parser_ast_node_positions) {
    auto p = parse("\n  int x = 1");
    const auto* d = as<VarDecl>(p.prog.stmts[0].get());
    CHECK(d && d->line == 2 && d->col == 3);
}

// ======================= TYPE CHECKER =======================
TEST(tc_valid_pipeline) {
    auto e = typeErrors("dna d = \"ATGAAA\"\nrna r = transcribe(d)\nprotein p = translate(r)\nprint(p)");
    CHECK(e.empty());
}

TEST(tc_dna_literal_invalid_char) {
    auto e = typeErrors("dna d = \"ACGX\"");
    CHECK_EQ(e.size(), size_t(1));
    CHECK(anyHas(e, "invalid character 'X'"));
}

TEST(tc_dna_literal_rejects_u) {
    CHECK(anyHas(typeErrors("dna d = \"ACGU\""), "invalid character 'U'"));
}

TEST(tc_rna_literal_rejects_t) {
    CHECK(anyHas(typeErrors("rna r = \"ACGT\""), "invalid character 'T'"));
}

TEST(tc_rna_literal_accepts_u) {
    CHECK(typeErrors("rna r = \"ACGU\"").empty());
}

TEST(tc_protein_alphabet) {
    CHECK(typeErrors("protein p = \"MKV*\"").empty());
    CHECK(anyHas(typeErrors("protein p = \"MKB\""), "invalid character 'B'"));  // B is not an amino acid code here
}

TEST(tc_empty_sequence_literal_ok) {
    CHECK(typeErrors("dna d = \"\"").empty());
}

TEST(tc_translate_requires_rna_with_hint) {
    auto e = typeErrors("dna d = \"ATG\"\nprotein p = translate(d)");
    CHECK_EQ(e.size(), size_t(1));
    CHECK(anyHas(e, "translate() expects rna, got dna"));
    CHECK(anyHas(e, "transcribe"));
}

TEST(tc_transcribe_requires_dna) {
    auto e = typeErrors("rna r = \"AUG\"\nrna x = transcribe(r)");
    CHECK(anyHas(e, "transcribe() expects dna, got rna"));
}

TEST(tc_revcomp_returns_dna) {
    CHECK(typeErrors("dna d = \"ATG\"\ndna r = revcomp(d)").empty());
    CHECK(!typeErrors("dna d = \"ATG\"\nrna r = revcomp(d)").empty());
}

TEST(tc_gc_content_accepts_dna_and_rna_returns_float) {
    CHECK(typeErrors("dna d = \"ATG\"\nfloat g = gc_content(d)").empty());
    CHECK(typeErrors("rna r = \"AUG\"\nfloat g = gc_content(r)").empty());
    CHECK(!typeErrors("dna d = \"ATG\"\nint g = gc_content(d)").empty());
}

TEST(tc_len_accepts_sequences_and_strings) {
    CHECK(typeErrors("dna d = \"AT\"\nint n = len(d)").empty());
    CHECK(typeErrors("string s = \"hello\"\nint n = len(s)").empty());
    CHECK(anyHas(typeErrors("int n = len(5)"), "len() expects"));
}

TEST(tc_string_literal_passed_to_function_rejected) {
    auto e = typeErrors("rna r = transcribe(\"ATG\")");
    CHECK(anyHas(e, "assign it to a typed variable first"));
}

TEST(tc_load_fasta_returns_dna) {
    CHECK(typeErrors("dna g = load_fasta(\"x.fa\")").empty());
    CHECK(typeErrors("string f = \"x.fa\"\ndna g = load_fasta(f)").empty());
    CHECK(!typeErrors("rna g = load_fasta(\"x.fa\")").empty());
}

TEST(tc_undefined_variable) {
    auto e = typeErrors("print(x)");
    CHECK_EQ(e.size(), size_t(1));
    CHECK(anyHas(e, "undefined variable 'x'"));
}

TEST(tc_no_cascade_after_undefined_variable) {
    auto e = typeErrors("rna r = transcribe(nope)");
    CHECK_EQ(e.size(), size_t(1));
}

TEST(tc_redeclaration) {
    auto e = typeErrors("int a = 1\nint a = 2");
    CHECK(anyHas(e, "already declared"));
}

TEST(tc_assignment_to_undeclared) {
    CHECK(anyHas(typeErrors("x = 5"), "undeclared variable 'x'"));
}

TEST(tc_assignment_type_mismatch) {
    CHECK(anyHas(typeErrors("int a = 1\na = \"hello\""), "cannot assign string to int"));
}

TEST(tc_int_widens_to_float_but_not_back) {
    CHECK(typeErrors("float f = 3").empty());
    CHECK(anyHas(typeErrors("int i = 2.5"), "cannot assign float to int"));
}

TEST(tc_string_variable_cannot_become_dna) {
    auto e = typeErrors("string s = \"ACGT\"\ndna d = s");
    CHECK(anyHas(e, "cannot assign string variable to dna"));
}

TEST(tc_string_literal_into_assign_checks_alphabet) {
    CHECK(anyHas(typeErrors("dna d = \"ACGT\"\nd = \"ACGZ\""), "invalid character 'Z'"));
}

TEST(tc_dna_to_rna_assignment_rejected) {
    CHECK(anyHas(typeErrors("dna d = \"ATG\"\nrna r = d"), "cannot assign dna to rna"));
}

TEST(tc_if_condition_must_be_bool) {
    CHECK(anyHas(typeErrors("if 5 {\nprint(1)\n}"), "if condition must be bool"));
    CHECK(typeErrors("if 5 > 3 {\nprint(1)\n}").empty());
}

TEST(tc_block_scope_in_if) {
    auto e = typeErrors("if 1 == 1 {\nint inner = 1\n}\nprint(inner)");
    CHECK(anyHas(e, "undefined variable 'inner'"));
}

TEST(tc_for_requires_list) {
    CHECK(anyHas(typeErrors("for x in 5 {\nprint(x)\n}"), "for loop expects list<int>"));
}

TEST(tc_for_loop_variable_is_int_and_scoped) {
    CHECK(typeErrors("dna g = \"AAAA\"\nfor i in find_all(/A/ in g) {\nint y = i + 1\n}").empty());
    auto e = typeErrors("dna g = \"AAAA\"\nfor i in find_all(/A/ in g) {\n}\nprint(i)");
    CHECK(anyHas(e, "undefined variable 'i'"));
}

TEST(tc_for_loop_variable_cannot_shadow) {
    auto e = typeErrors("int i = 0\ndna g = \"AA\"\nfor i in find_all(/A/ in g) {\n}");
    CHECK(anyHas(e, "loop variable 'i' is already declared"));
}

TEST(tc_find_all_ok) {
    CHECK(typeErrors("dna g = \"ATAT\"\nmotif m = /AT/\nint n = count(find_all(m in g))").empty());
}

TEST(tc_find_all_works_on_rna_and_protein) {
    CHECK(typeErrors("rna r = \"AUAU\"\nint n = count(find_all(/AU/ in r))").empty());
    CHECK(typeErrors("protein p = \"MKV\"\nint n = count(find_all(/K/ in p))").empty());
}

TEST(tc_find_all_wrong_haystack) {
    CHECK(anyHas(typeErrors("int n = count(find_all(/A/ in 5))"), "find_all() expects dna, rna or protein after 'in'"));
}

TEST(tc_find_all_wrong_needle) {
    auto e = typeErrors("dna g = \"AT\"\nint n = count(find_all(g in g))");
    CHECK(anyHas(e, "find_all() expects a motif before 'in'"));
}

TEST(tc_find_all_missing_in) {
    CHECK(anyHas(typeErrors("dna g = \"AT\"\nint n = count(find_all(g))"), "motif in sequence"));
}

TEST(tc_in_outside_find_all) {
    CHECK(anyHas(typeErrors("dna g = \"AT\"\nprint(len(/A/ in g))"), "'in' can only be used inside find_all()"));
}

TEST(tc_match_end_ok) {
    CHECK(typeErrors("dna g = \"ATAT\"\nint e = match_end(/AT/ in g, 0)").empty());
}

TEST(tc_match_end_position_must_be_int) {
    CHECK(anyHas(typeErrors("dna g = \"AT\"\nint e = match_end(/A/ in g, 1.5)"), "position must be int"));
}

TEST(tc_match_end_wrong_arity) {
    CHECK(anyHas(typeErrors("dna g = \"AT\"\nint e = match_end(/A/ in g)"), "match_end() expects"));
}

TEST(tc_count_requires_list) {
    CHECK(anyHas(typeErrors("int n = count(5)"), "count() expects list<int>"));
}

TEST(tc_wrong_arity) {
    CHECK(anyHas(typeErrors("dna d = \"A\"\nint n = len(d, d)"), "len() takes 1 argument, got 2"));
    CHECK(anyHas(typeErrors("int n = len()"), "len() takes 1 argument, got 0"));
}

TEST(tc_unknown_function) {
    CHECK(anyHas(typeErrors("int n = frobnicate(1)"), "unknown function 'frobnicate'"));
}

TEST(tc_slice_keeps_type) {
    CHECK(typeErrors("dna d = \"ACGT\"\ndna s = d[0:2]").empty());
    CHECK(typeErrors("rna r = \"ACGU\"\nrna s = r[1:]").empty());
    CHECK(typeErrors("protein p = \"MKV\"\nprotein s = p[:2]").empty());
    CHECK(typeErrors("string t = \"hello\"\nstring s = t[1:3]").empty());
}

TEST(tc_slice_type_mismatch) {
    CHECK(anyHas(typeErrors("dna d = \"ACGT\"\nprotein s = d[0:2]"), "cannot assign dna to protein"));
}

TEST(tc_slice_bound_must_be_int) {
    CHECK(anyHas(typeErrors("dna d = \"ACGT\"\ndna s = d[0.5:2]"), "slice bound must be int"));
}

TEST(tc_cannot_slice_int) {
    CHECK(anyHas(typeErrors("int a = 5\nprint(a[0:1])"), "cannot slice int"));
}

TEST(tc_list_indexing) {
    CHECK(typeErrors("dna g = \"AAAA\"\nint x = find_all(/A/ in g)[0]").empty());
    CHECK(anyHas(typeErrors("dna g = \"AAAA\"\nint x = find_all(/A/ in g)[\"a\"]"), "index must be int"));
}

TEST(tc_cannot_index_sequence_directly_with_hint) {
    auto e = typeErrors("dna d = \"ACGT\"\nprint(d[0])");
    CHECK(anyHas(e, "cannot index dna"));
    CHECK(anyHas(e, "use a slice"));
}

TEST(tc_cannot_index_int) {
    auto e = typeErrors("int a = 5\nprint(a[0])");
    CHECK(anyHas(e, "cannot index int"));
    CHECK(!anyHas(e, "use a slice"));
}

TEST(tc_concatenation) {
    CHECK(typeErrors("dna a = \"AC\"\ndna b = \"GT\"\ndna c = a + b").empty());
    CHECK(typeErrors("string a = \"x\"\nstring b = \"y\"\nstring c = a + b").empty());
    CHECK(typeErrors("rna a = \"AC\"\nrna b = \"GU\"\nrna c = a + b").empty());
}

TEST(tc_concatenation_mismatch) {
    auto e = typeErrors("dna a = \"AC\"\nrna b = \"GU\"\nprint(a + b)");
    CHECK(anyHas(e, "operator '+' cannot be applied to dna and rna"));
}

TEST(tc_sequence_subtraction_rejected) {
    CHECK(anyHas(typeErrors("dna a = \"AC\"\ndna b = \"GT\"\nprint(a - b)"), "operator '-' cannot be applied"));
}

TEST(tc_arithmetic_result_types) {
    CHECK(typeErrors("int a = 1 + 2 * 3").empty());
    CHECK(typeErrors("float a = 1 + 2.5").empty());
    CHECK(!typeErrors("int a = 1 + 2.5").empty());
    CHECK(typeErrors("float a = 7 / 2").empty());     // int / int is int, widened to float
}

TEST(tc_equality_rules) {
    CHECK(typeErrors("bool b = 1 == 1.0").empty());
    CHECK(typeErrors("dna a = \"A\"\ndna b = \"A\"\nbool c = a == b").empty());
    auto e = typeErrors("dna a = \"A\"\nrna b = \"A\"\nbool c = a == b");
    CHECK(anyHas(e, "cannot compare dna with rna"));
}

TEST(tc_ordering_only_for_numbers) {
    CHECK(typeErrors("bool b = 1 < 2.5").empty());
    CHECK(anyHas(typeErrors("dna a = \"A\"\ndna b = \"C\"\nbool c = a < b"), "operator '<' cannot be applied"));
}

TEST(tc_bool_from_comparison_in_if) {
    CHECK(typeErrors("dna d = \"ATG\"\nif len(d) >= 3 {\nprint(1)\n}").empty());
}

TEST(tc_invalid_motif_reports_error_with_column) {
    auto e = typeErrors("motif m = /(AT/");
    CHECK_EQ(e.size(), size_t(1));
    CHECK(anyHas(e, "invalid motif"));
    CHECK(e[0].find("error at 1:") == 0);
}

TEST(tc_invalid_motif_variants) {
    for (const char* pat : {"[AC", "A{2", "A{3,1}", "*A", "(A))"}) {
        auto e = typeErrors(std::string("motif m = /") + pat + "/");
        CHECK(anyHas(e, "invalid motif"));
    }
}

TEST(tc_valid_motif_syntax_all_accepted) {
    for (const char* pat : {"ATG", "A.C", "[ACG]T", "[^AT]", "[A-C]", "A|C", "A|", "(AT)*", "A+", "A?",
                            "A{3}", "A{2,}", "A{2,4}", "A+?", "ATG(<sense>){3,}<stop>", "<codon>"}) {
        auto e = typeErrors(std::string("motif m = /") + pat + "/");
        if (!e.empty()) std::printf("    (pattern /%s/ gave: %s)\n", pat, e[0].c_str());
        CHECK(e.empty());
    }
}

TEST(tc_motifs_are_collected) {
    Parsed p = parse("motif a = /AT/\nmotif b = /GC+/\ndna g = \"ATGC\"\nint n = count(find_all(/T/ in g))");
    CHECK(p.errors.empty());
    TypeChecker tc;
    auto e = tc.check(p.prog);
    CHECK(e.empty());
    CHECK_EQ(tc.motifs().size(), size_t(3));
    for (const auto& m : tc.motifs()) CHECK(m != nullptr);
}

TEST(tc_motif_dfa_attached_to_ast) {
    Parsed p = parse("motif m = /AT+/");
    TypeChecker tc;
    tc.check(p.prog);
    const auto* d = as<VarDecl>(p.prog.stmts[0].get());
    const auto* r = as<RegexLit>(d->init.get());
    CHECK(r && r->dfa != nullptr);
    if (r && r->dfa) CHECK(r->dfa->numStates() > 0);
}

TEST(tc_checker_is_reusable) {
    TypeChecker tc;
    Parsed bad = parse("print(x)");
    CHECK(!tc.check(bad.prog).empty());
    Parsed good = parse("int x = 1\nprint(x)");
    CHECK(tc.check(good.prog).empty());     // state from the first run must not leak
}

TEST(tc_multiple_errors_all_reported) {
    auto e = typeErrors("dna a = \"AX\"\nrna b = \"AT\"\nint c = \"s\"");
    CHECK_EQ(e.size(), size_t(3));
}

TEST(tc_error_positions_are_one_based_and_correct) {
    auto e = typeErrors("int a = 1\nprint(zzz)");
    CHECK_EQ(e.size(), size_t(1));
    CHECK(e[0].find("error at 2:7") == 0);
}

// ---------- runner ----------
int main(int argc, char** argv) {
    std::string filter = argc > 1 ? argv[1] : "";
    int ran = 0, failedTests = 0;
    for (const auto& t : registry()) {
        if (!filter.empty() && std::string(t.name).find(filter) == std::string::npos) continue;
        g_current = t.name;
        int before = g_failedChecks;
        t.fn();
        ran++;
        if (g_failedChecks != before) { failedTests++; std::printf("FAIL  %s\n", t.name); }
    }
    std::printf("\n%d tests run, %d passed, %d failed\n", ran, ran - failedTests, failedTests);
    return failedTests == 0 ? 0 : 1;
}
