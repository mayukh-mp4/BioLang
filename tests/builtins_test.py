"""Tests for the interpreter built-ins: transcribe, translate, revcomp, gc_content, len,
count, load_fasta, concatenation, comparisons, control flow, runtime errors.

Usage (from the project root):
    python tests\\builtins_test.py build\\biolang.exe

The reference genetic code below is written out by amino acid (independent of the table in
interpreter.cpp, which is ordered T,C,A,G). Exhaustive checks cover all 64 codons and 300 random
sequences for revcomp / gc_content / translate / transcribe.
"""
import itertools, os, random, subprocess, sys, tempfile

exe = sys.argv[1] if len(sys.argv) > 1 else ("build/biolang.exe" if os.name == "nt" else "build/biolang")
tmp = tempfile.mkdtemp()

AA = {  # amino acid -> codons (standard code, DNA)
    "F": "TTT TTC", "L": "TTA TTG CTT CTC CTA CTG", "S": "TCT TCC TCA TCG AGT AGC",
    "Y": "TAT TAC", "*": "TAA TAG TGA", "C": "TGT TGC", "W": "TGG",
    "P": "CCT CCC CCA CCG", "H": "CAT CAC", "Q": "CAA CAG",
    "R": "CGT CGC CGA CGG AGA AGG", "I": "ATT ATC ATA", "M": "ATG",
    "T": "ACT ACC ACA ACG", "N": "AAT AAC", "K": "AAA AAG",
    "V": "GTT GTC GTA GTG", "A": "GCT GCC GCA GCG", "D": "GAT GAC", "E": "GAA GAG",
    "G": "GGT GGC GGA GGG",
}
CODON = {c: a for a, cs in AA.items() for c in cs.split()}
assert len(CODON) == 64

def ref_translate(dna):
    return "".join(CODON[dna[i:i + 3]] for i in range(0, len(dna) - 2, 3))

def ref_revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))

def ref_gc(s):
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0

passed = failed = 0

def run(src, extra=()):
    p = os.path.join(tmp, "t.bl")
    with open(p, "w", newline="\n") as f:
        f.write(src)
    r = subprocess.run([exe, *extra, p], capture_output=True, text=True)
    return r.returncode, r.stdout.strip(), r.stderr.strip()

def check(name, src, want, rc=0, err=None):
    global passed, failed
    code, out, e = run(src)
    ok = code == rc and (out == want if want is not None else True) and (err is None or err in e)
    if ok: passed += 1
    else:
        failed += 1
        print(f"FAIL  {name}\n   want rc={rc} out={want!r} err~{err!r}\n   got  rc={code} out={out!r} err={e[:150]!r}")

# ---- directed cases ----
check("transcribe", 'dna d = "ATGCCT"\nprint(transcribe(d))', "AUGCCU")
check("translate", 'dna d = "ATGGCCTAA"\nprint(translate(transcribe(d)))', "MA*")
check("translate drops partial codon", 'dna d = "ATGGC"\nprint(translate(transcribe(d)))', "M")
check("translate empty", 'dna d = ""\nprint(translate(transcribe(d)))', "")
check("revcomp", 'dna d = "AACG"\nprint(revcomp(d))', "CGTT")
check("revcomp palindrome", 'dna d = "GAATTC"\nprint(revcomp(d) == d)', "true")
check("gc_content", 'dna d = "GGCCAATT"\nprint(gc_content(d))', "0.5")
check("gc_content empty", 'dna d = ""\nprint(gc_content(d))', "0")
check("gc_content rna", 'rna r = "GCAU"\nprint(gc_content(r))', "0.5")
check("len", 'dna d = "ACGTA"\nprint(len(d))', "5")
check("len string", 'string s = "hello"\nprint(len(s))', "5")
check("count", 'dna d = "AAAA"\nprint(count(find_all(/A/ in d)))', "4")
check("concat dna", 'dna a = "AC"\ndna b = "GT"\nprint(a + b)', "ACGT")
check("concat then len", 'dna a = "AC"\ndna b = "GT"\nprint(len(a + b))', "4")
check("slice", 'dna d = "ACGTAC"\nprint(d[1:4])', "CGT")
check("slice open ends", 'dna d = "ACGT"\nprint(d[:2])\nprint(d[2:])\nprint(d[:])', "AC\nGT\nACGT")
check("slice empty", 'dna d = "ACGT"\nprint(len(d[2:2]))', "0")
check("slice out of range", 'dna d = "ACGT"\nprint(d[0:9])', None, rc=1, err="out of range")
check("slice reversed", 'dna d = "ACGT"\nprint(d[3:1])', None, rc=1, err="out of range")
check("list index", 'dna d = "ACAC"\nprint(find_all(/A/ in d)[1])', "2")
check("list index out of range", 'dna d = "ACAC"\nprint(find_all(/A/ in d)[5])', None, rc=1, err="out of range")
check("match_end hit", 'dna d = "GGATGAA"\nprint(match_end(/ATG/ in d, 2))', "5")
check("match_end miss", 'dna d = "GGATGAA"\nprint(match_end(/ATG/ in d, 0))', "-1")
check("match_end bad pos", 'dna d = "GG"\nprint(match_end(/G/ in d, 9))', None, rc=1, err="out of range")
check("division by zero", "print(1 / 0)", None, rc=1, err="division by zero")
check("float division by zero", "print(1.0 / 0.0)", None, rc=1, err="division by zero")
check("int arithmetic", "print(7 / 2)\nprint(7 - 2 * 3)\nprint(2 + 3 * 4)", "3\n1\n14")
check("float arithmetic", "print(7 / 2.0)", "3.5")
check("comparisons", "print(1 < 2)\nprint(2 <= 2)\nprint(3 > 4)\nprint(2 >= 3)\nprint(1 == 1.0)\nprint(1 != 2)",
      "true\ntrue\nfalse\nfalse\ntrue\ntrue")
check("sequence equality", 'dna a = "AC"\ndna b = "AC"\ndna c = "AG"\nprint(a == b)\nprint(a == c)', "true\nfalse")
check("if else", "if 1 > 2 {\nprint(1)\n} else {\nprint(2)\n}", "2")
check("else if", "int x = 5\nif x < 3 {\nprint(1)\n} else if x < 10 {\nprint(2)\n} else {\nprint(3)\n}", "2")
check("for loop", 'dna d = "AAAA"\nfor i in find_all(/A/ in d) {\nprint(i)\n}', "0\n1\n2\n3")
check("for empty list", 'dna d = "CCCC"\nfor i in find_all(/A/ in d) {\nprint(i)\n}\nprint(9)', "9")
check("nested loop+if", 'dna d = "ACACAC"\nint s = 0\nfor i in find_all(/A/ in d) {\nif i > 0 {\ns = s + i\n}\n}\nprint(s)', "6")
check("assignment", "int x = 1\nx = x + 4\nprint(x)", "5")
check("int widens to float", "float f = 3\nprint(f / 2)", "1.5")
check("list printing", 'dna d = "ACAC"\nprint(find_all(/A/ in d))', "[0, 2]")
check("empty list printing", 'dna d = "CC"\nprint(find_all(/A/ in d))', "[]")
check("motif printing", "motif m = /AT+/\nprint(m)", "/AT+/")
check("bool printing", "print(1 == 1)", "true")
check("type error blocks execution", 'dna d = "AX"\nprint(1)', None, rc=1, err="invalid character")
check("lex error", "x = @", None, rc=1, err="lexical error")
check("parse error", "int x =", None, rc=1, err="expected an expression")

# ---- load_fasta ----
def fasta(name, text):
    p = os.path.join(tmp, name).replace("\\", "/")
    with open(p, "w", newline="\n") as f: f.write(text)
    return p

f1 = fasta("a.fa", ">rec1 description\nacgt\nACGT\n\n>rec2\nTTTT\n")
check("fasta: first record, uppercased, joined", f'dna g = load_fasta("{f1}")\nprint(g)\nprint(len(g))', "ACGTACGT\n8")
f2 = fasta("b.fa", ">x\r\nACGT\r\nNNRY\r\n")
check("fasta: CRLF and IUPAC codes", f'dna g = load_fasta("{f2}")\nprint(g)', "ACGTNNRY")
f3 = fasta("c.fa", ">x\nACGJ\n")
check("fasta: invalid character", f'dna g = load_fasta("{f3}")', None, rc=1, err="invalid character 'J'")
f4 = fasta("d.fa", ">x\n")
check("fasta: no sequence", f'dna g = load_fasta("{f4}")', None, rc=1, err="no sequence found")
check("fasta: missing file", 'dna g = load_fasta("/nonexistent/zzz.fa")', None, rc=1, err="cannot open FASTA")
check("fasta: path in string variable", f'string p = "{f1}"\ndna g = load_fasta(p)\nprint(len(g))', "8")
f5 = fasta("e.fa", ">x\nACGT\nACGT")   # no trailing newline
check("fasta: no trailing newline", f'dna g = load_fasta("{f5}")\nprint(len(g))', "8")

# ---- exhaustive / randomized vs reference ----
bad = 0
for c in map("".join, itertools.product("ACGT", repeat=3)):
    _, out, _ = run(f'dna d = "{c}"\nprint(translate(transcribe(d)))')
    if out != CODON[c]:
        bad += 1; print(f"FAIL  codon {c}: got {out!r} want {CODON[c]!r}")
passed += 64 - bad; failed += bad

random.seed(7)
bad = 0
N = 300
for _ in range(N):
    n = random.randint(1, 60)
    s = "".join(random.choice("ACGT") for _ in range(n))
    prog = (f'dna d = "{s}"\nprint(revcomp(d))\nprint(transcribe(d))\nprint(translate(transcribe(d)))\n'
            f'print(gc_content(d))\nprint(len(d))\nprint(revcomp(revcomp(d)) == d)')
    _, out, _ = run(prog)
    lines = out.split("\n")
    want = [ref_revcomp(s), s.replace("T", "U"), ref_translate(s), None, str(n), "true"]
    ok = len(lines) == 6 and all(w is None or l == w for l, w in zip(lines, want))
    if ok:
        ok = abs(float(lines[3]) - ref_gc(s)) < 1e-5     # cout prints 6 significant digits
    if not ok:
        bad += 1; print(f"FAIL  random seq {s}\n   got {lines}")
passed += N - bad; failed += bad

print(f"\n{passed} passed, {failed} failed")
print("ALL PASS" if not failed else "SOME FAILED")
sys.exit(1 if failed else 0)
