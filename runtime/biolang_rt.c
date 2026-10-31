/* BioLang runtime for the LLVM backend.
 * Link with the .ll file produced by `biolang --emit-llvm`:
 *     clang prog.ll runtime/biolang_rt.c -O2 -o prog
 * Mirrors the interpreter's semantics and error messages. Memory is never freed
 * (programs are short-lived scripts). */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { long long n; char d[]; } Str;            /* layout shared with the generated IR */
typedef struct { long long n, cap; long long *d; } List;
typedef struct Motif {
    const Str *pattern;
    int K, start, shortest;
    const unsigned char *classMap;
    const int *trans;
    const unsigned char *acc;
    const unsigned char *startByte;
    const int *rtrans;                                    /* NULL if the reverse DFA was too large */
    const unsigned char *racc;
    long long (*match)(const char *, long long, long long, long long *);
} Motif;

static void fail(int line, int col, const char *fmt, ...) __attribute__((noreturn, format(printf, 3, 4)));
#include <stdarg.h>
static void fail(int line, int col, const char *fmt, ...) {
    char msg[512];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(msg, sizeof msg, fmt, ap);
    va_end(ap);
    fflush(stdout);
    fprintf(stderr, "runtime error at %d:%d: %s\n", line, col, msg);
    exit(1);
}

static Str *new_str(long long n) {
    Str *s = (Str *)malloc(sizeof(Str) + (size_t)n + 1);
    if (!s) { fprintf(stderr, "out of memory\n"); exit(1); }
    s->n = n;
    s->d[n] = 0;
    return s;
}

/* ---------- arithmetic checks ---------- */
void bl_chk_div_i(long long d, int line, int col) { if (d == 0) fail(line, col, "division by zero"); }
void bl_chk_div_f(double d, int line, int col)    { if (d == 0.0) fail(line, col, "division by zero"); }

/* ---------- strings / sequences ---------- */
long long bl_len(const Str *s) { return s->n; }

Str *bl_concat(const Str *a, const Str *b) {
    Str *r = new_str(a->n + b->n);
    memcpy(r->d, a->d, (size_t)a->n);
    memcpy(r->d + a->n, b->d, (size_t)b->n);
    return r;
}

Str *bl_slice(const Str *s, long long lo, long long hi, int line, int col) {
    if (lo < 0 || hi < 0 || lo > s->n || hi > s->n || lo > hi)
        fail(line, col, "slice [%lld:%lld] out of range (length %lld)", lo, hi, s->n);
    Str *r = new_str(hi - lo);
    memcpy(r->d, s->d + lo, (size_t)(hi - lo));
    return r;
}

int bl_str_eq(const Str *a, const Str *b) { return a->n == b->n && memcmp(a->d, b->d, (size_t)a->n) == 0; }

static const char *CODON_TABLE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

static int base_index(char c) {
    switch (c) { case 'T': case 'U': return 0; case 'C': return 1; case 'A': return 2; case 'G': return 3; }
    return -1;
}

Str *bl_transcribe(const Str *s) {
    Str *r = new_str(s->n);
    for (long long i = 0; i < s->n; i++) r->d[i] = s->d[i] == 'T' ? 'U' : s->d[i];
    return r;
}

Str *bl_translate(const Str *s) {
    Str *r = new_str(s->n / 3);
    long long k = 0;
    for (long long i = 0; i + 2 < s->n; i += 3) {
        int a = base_index(s->d[i]), b = base_index(s->d[i + 1]), c = base_index(s->d[i + 2]);
        r->d[k++] = (a < 0 || b < 0 || c < 0) ? 'X' : CODON_TABLE[16 * a + 4 * b + c];
    }
    r->n = k;
    r->d[k] = 0;
    return r;
}

Str *bl_revcomp(const Str *s) {
    Str *r = new_str(s->n);
    for (long long i = 0; i < s->n; i++) {
        char c = s->d[s->n - 1 - i];
        switch (c) {
            case 'A': c = 'T'; break; case 'T': c = 'A'; break;
            case 'C': c = 'G'; break; case 'G': c = 'C'; break;
            case 'R': c = 'Y'; break; case 'Y': c = 'R'; break;
            case 'K': c = 'M'; break; case 'M': c = 'K'; break;
            case 'B': c = 'V'; break; case 'V': c = 'B'; break;
            case 'D': c = 'H'; break; case 'H': c = 'D'; break;
        }
        r->d[i] = c;
    }
    return r;
}

double bl_gc(const Str *s) {
    if (s->n == 0) return 0.0;
    long long gc = 0;
    for (long long i = 0; i < s->n; i++) if (s->d[i] == 'G' || s->d[i] == 'C') gc++;
    return (double)gc / (double)s->n;
}

/* First record of a FASTA file; same rules as the interpreter. */
Str *bl_load_fasta(const Str *path, int line, int col) {
    FILE *f = fopen(path->d, "rb");
    if (!f) fail(line, col, "cannot open FASTA file '%s'", path->d);
    static const char *allowed = "ACGTNRYKMSWBDHV";
    size_t cap = 1 << 16, n = 0;
    char *seq = (char *)malloc(cap);
    int inRecord = 0, stop = 0;
    char *ln = NULL;
    size_t lcap = 0;
    int c;
    while (!stop) {
        size_t len = 0;
        while ((c = fgetc(f)) != EOF && c != '\n') {
            if (len + 1 >= lcap) { lcap = lcap ? lcap * 2 : 1024; ln = (char *)realloc(ln, lcap); }
            ln[len++] = (char)c;
        }
        if (c == EOF && len == 0) break;
        if (len && ln[len - 1] == '\r') len--;
        if (len == 0 || ln[0] == ';') { if (c == EOF) break; continue; }
        if (ln[0] == '>') {
            if (inRecord) break;
            inRecord = 1;
            if (c == EOF) break;
            continue;
        }
        inRecord = 1;
        for (size_t i = 0; i < len; i++) {
            if (ln[i] == ' ') continue;
            char u = (char)toupper((unsigned char)ln[i]);
            if (!strchr(allowed, u) || u == 0) {
                fclose(f);
                fail(line, col, "invalid character '%c' in FASTA file '%s'", ln[i], path->d);
            }
            if (n + 1 >= cap) { cap *= 2; seq = (char *)realloc(seq, cap); }
            seq[n++] = u;
        }
        if (c == EOF) break;
    }
    fclose(f);
    if (n == 0) fail(line, col, "no sequence found in '%s'", path->d);
    Str *r = new_str((long long)n);
    memcpy(r->d, seq, n);
    free(seq);
    free(ln);
    return r;
}

/* ---------- lists ---------- */
static List *new_list(void) { return (List *)calloc(1, sizeof(List)); }
static void push(List *l, long long v) {
    if (l->n == l->cap) { l->cap = l->cap ? l->cap * 2 : 16; l->d = (long long *)realloc(l->d, (size_t)l->cap * sizeof(long long)); }
    l->d[l->n++] = v;
}
long long bl_list_len(const List *l) { return l->n; }
long long bl_list_at(const List *l, long long i) { return l->d[i]; }
long long bl_list_get(const List *l, long long i, int line, int col) {
    if (i < 0 || i >= l->n) fail(line, col, "index %lld out of range (list has %lld items)", i, l->n);
    return l->d[i];
}

/* ---------- motif scanning (same strategy as Dfa::findAll in Auto mode) ---------- */
static int scan_restart(const Motif *m, const char *s, long long n, List *out, long long budget) {
    long long work = 0, i = 0;
    while (i < n) {
        while (i < n && !m->startByte[(unsigned char)s[i]]) i++;
        if (i >= n) break;
        long long end = m->match(s, n, i, &work);
        if (budget >= 0 && work > budget) return 1;
        if (end > i) { push(out, i); i = end; } else i++;
    }
    return 0;
}

static void scan_linear(const Motif *m, const char *s, long long n, List *out) {
    unsigned char *isStart = (unsigned char *)malloc((size_t)n);
    int state = 0;
    for (long long i = n; i-- > 0;) {
        state = m->rtrans[state * m->K + m->classMap[(unsigned char)s[i]]];
        isStart[i] = m->racc[state];
    }
    long long work = 0, i = 0;
    while (i < n) {
        while (i < n && !isStart[i]) i++;
        if (i >= n) break;
        long long end = m->match(s, n, i, &work);
        if (end > i) { push(out, i); i = end; } else i++;
    }
    free(isStart);
}

List *bl_find_all(const Motif *m, const Str *seq) {
    List *out = new_list();
    long long n = seq->n;
    if (n == 0) return out;
    if (scan_restart(m, seq->d, n, out, 16 * n + 4096) == 0) return out;
    out->n = 0;                                           /* budget exceeded: adversarial input */
    if (m->rtrans) { scan_linear(m, seq->d, n, out); return out; }
    scan_restart(m, seq->d, n, out, -1);
    return out;
}

long long bl_match_end(const Motif *m, const Str *seq, long long p, int line, int col) {
    if (p < 0 || p > seq->n) fail(line, col, "match_end position %lld out of range", p);
    long long work = 0;
    return m->match(seq->d, seq->n, p, &work);
}

/* ---------- printing ---------- */
void bl_print_int(long long v)   { printf("%lld\n", v); }
void bl_print_float(double v)    { printf("%g\n", v); }
void bl_print_bool(int v)        { puts(v ? "true" : "false"); }
void bl_print_str(const Str *s)  { fwrite(s->d, 1, (size_t)s->n, stdout); putchar('\n'); }
void bl_print_motif(const Motif *m) { printf("/%s/\n", m->pattern->d); }
void bl_print_list(const List *l) {
    putchar('[');
    for (long long i = 0; i < l->n; i++) printf("%s%lld", i ? ", " : "", l->d[i]);
    puts("]");
}
