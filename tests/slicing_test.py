"""Tests for slicing, string variables, list indexing, concatenation and match_end.
usage: python tests/slicing_test.py build/biolang.exe"""
import subprocess, sys, tempfile, os, random

exe = sys.argv[1]
fails = 0

def run(src):
    with tempfile.NamedTemporaryFile("w", suffix=".bl", delete=False) as f:
        f.write(src)
    p = subprocess.run([exe, f.name], capture_output=True, text=True)
    os.unlink(f.name)
    return p.returncode, p.stdout.strip(), p.stderr.strip()

def ok(name, src, expect_out):
    global fails
    rc, out, err = run(src)
    good = rc == 0 and out == expect_out
    print(("PASS " if good else "FAIL ") + name)
    if not good:
        fails += 1; print("   got:", rc, repr(out), err)

def bad(name, src, needle, runtime=False):
    global fails
    rc, out, err = run(src)
    good = rc != 0 and needle in err
    print(("PASS " if good else "FAIL ") + name)
    if not good:
        fails += 1; print("   got:", rc, repr(out), err)

ok("slice basic",      'dna d = "ACGTACGTTT"\nprint(d[2:6])', "GTAC")
ok("slice open ends",  'dna d = "ACGTACGTTT"\nprint(d[:2])\nprint(d[8:])\nprint(d[:])', "AC\nTT\nACGTACGTTT")
ok("empty slice",      'dna d = "ACGT"\nprint(len(d[2:2]))', "0")
ok("slice keeps type", 'dna d = "ACGTACGTA"\nprint(translate(transcribe(d[0:6])))', "TY")
ok("expr bounds",      'dna d = "ACGTACGTTT"\nint k = 3\nprint(d[k-1:k*2])', "GTAC")
ok("chained slice",    'dna d = "ACGTACGTTT"\nprint(d[2:8][1:3])', "TA")
ok("concat",           'dna a = "AC"\ndna b = "GT"\nprint(a + b + a)', "ACGTAC")
ok("string var",       'string s = "hello"\nprint(s[1:3])\nprint(len(s))', "el\n5")
ok("division after ]", 'dna d = "ACGTACGT"\nprint(len(d[0:4]) / 2)', "2")
ok("list index",       'dna d = "ACGTACGT"\nmotif m = /ACG/\nprint(find_all(m in d)[1])', "4")
ok("match_end",        'dna d = "TTATGAAACCCTAAT"\nmotif m = /ATG(<sense>)*<stop>/\nprint(match_end(m in d, 2))', "14")
ok("match_end miss",   'dna d = "ACGT"\nmotif m = /TT/\nprint(match_end(m in d, 0))', "-1")
open("_t.fa", "w").write(">x\nACGTACGTAC\n")
ok("fasta via string var", 'string f = "_t.fa"\ndna g = load_fasta(f)\nprint(g[2:6])', "GTAC")
os.unlink("_t.fa")

bad("slice out of range", 'dna d = "ACGT"\nprint(d[0:9])', "out of range")
bad("slice start>end",    'dna d = "ACGT"\nprint(d[3:1])', "out of range")
bad("negative bound",     'dna d = "ACGT"\nprint(d[0-1:2])', "out of range")
bad("list index range",   'dna d = "ACGT"\nmotif m = /ACG/\nprint(find_all(m in d)[5])', "out of range")
bad("float bound",        'dna d = "ACGT"\nprint(d[0:1.5])', "slice bound must be int")
bad("slice an int",       'int x = 5\nprint(x[0:1])', "cannot slice int")
bad("index a dna",        'dna d = "ACGT"\nprint(d[1])', "use a slice")
bad("string var -> dna",  'string s = "ACGT"\ndna d = s', "cannot assign string variable")
bad("mix dna+rna",        'dna d = "ACGT"\nrna r = "ACGU"\nprint(d + r)', "cannot be applied")
bad("bad literal in string var kept", 'dna d = "ACXT"', "invalid character")
bad("unclosed bracket",   'dna d = "ACGT"\nprint(d[1:2)', "expected ']'")
bad("match_end wrong form", 'dna d = "ACGT"\nprint(match_end(d, 1))', "match_end() expects")

# fuzz: random slices must agree with Python slicing
random.seed(7)
mism = 0
for _ in range(300):
    n = random.randint(0, 30)
    s = "".join(random.choice("ACGT") for _ in range(n))
    lo = random.randint(0, n); hi = random.randint(lo, n)
    rc, out, err = run(f'dna d = "{s}"\nprint(d[{lo}:{hi}])' if n else 'print(0)')
    if n and out != s[lo:hi]: mism += 1
print("fuzz mismatches:", mism)
fails += mism
print("ALL PASS" if fails == 0 else f"{fails} FAILURE(S)")
sys.exit(1 if fails else 0)
