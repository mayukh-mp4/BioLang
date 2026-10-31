"""Differential test: the LLVM backend must behave exactly like the interpreter.

Every program is run twice, by `biolang prog.bl` (tree-walking interpreter) and as a native binary
(`biolang --emit-llvm` -> clang -> run). stdout, stderr and the exit code must be identical.

Usage (from the project root):
    python tests\\llvm_diff_test.py build\\biolang.exe

Needs clang (LLVM >= 15) on PATH; set CLANG to use another binary.
Program sources come from the directed cases in tests/builtins_test.py plus extra cases below:
motif patterns on 200 kb random sequences, adversarial inputs (budget -> linear fallback), and
random-pattern fuzzing.
"""
import ast, os, random, subprocess, sys, tempfile

exe = sys.argv[1] if len(sys.argv) > 1 else ("build/biolang.exe" if os.name == "nt" else "build/biolang")
here = os.path.dirname(os.path.abspath(__file__))
rt_src = os.path.join(here, "..", "runtime", "biolang_rt.c")
clang = os.environ.get("CLANG", "clang")
tmp = tempfile.mkdtemp()
CFLAGS = ["-O1", "-Wno-override-module"]
NATIVE_TIMEOUT = 5          # seconds; the interpreter needs ~10 ms on the adversarial cases

rt_obj = os.path.join(tmp, "rt.o")
subprocess.run([clang, "-O2", "-c", rt_src, "-o", rt_obj], check=True)

passed = failed = 0
bin_ext = ".exe" if os.name == "nt" else ""
counter = 0


def both(src):
    """Returns (interpreter result, native result) as (rc, stdout, stderr) tuples."""
    global counter
    counter += 1
    bl = os.path.join(tmp, f"p{counter}.bl")
    with open(bl, "w", newline="\n") as f:
        f.write(src)
    i = subprocess.run([exe, bl], capture_output=True, text=True)
    interp = (i.returncode, i.stdout, i.stderr.replace(bl, "<f>"))
    e = subprocess.run([exe, "--emit-llvm", bl], capture_output=True, text=True)
    if e.returncode != 0:                                  # front-end error: same message expected
        return interp, (e.returncode, "", e.stderr.replace(bl, "<f>"))
    ll = os.path.join(tmp, f"p{counter}.ll")
    with open(ll, "w", newline="\n") as f:
        f.write(e.stdout)
    out = os.path.join(tmp, f"p{counter}{bin_ext}")
    c = subprocess.run([clang, *CFLAGS, ll, rt_obj, "-o", out], capture_output=True, text=True)
    if c.returncode != 0:
        return interp, (-999, "", "CLANG FAILED: " + c.stderr[:400])
    try:
        r = subprocess.run([out], capture_output=True, text=True, timeout=NATIVE_TIMEOUT)
    except subprocess.TimeoutExpired:
        return interp, (-998, "", f"NATIVE TIMED OUT (> {NATIVE_TIMEOUT}s): linear-scan fallback lost?")
    return interp, (r.returncode, r.stdout, r.stderr)


def check(name, src):
    global passed, failed
    a, b = both(src)
    if a == b:
        passed += 1
    else:
        failed += 1
        print(f"FAIL  {name}")
        print(f"   interpreter: rc={a[0]} out={a[1][:150]!r} err={a[2][:150]!r}")
        print(f"   native     : rc={b[0]} out={b[1][:150]!r} err={b[2][:150]!r}")


# ---- 1. directed programs, harvested from builtins_test.py (string-constant cases only) ----
tree = ast.parse(open(os.path.join(here, "builtins_test.py")).read())
directed = []
for node in ast.walk(tree):
    if (isinstance(node, ast.Call) and getattr(node.func, "id", "") == "check" and len(node.args) >= 2
            and isinstance(node.args[0], ast.Constant) and isinstance(node.args[1], ast.Constant)):
        directed.append((node.args[0].value, node.args[1].value))
for name, src in directed:
    check("directed: " + name, src)

# ---- 2. extra directed programs (FASTA, scoping, mixed arithmetic, nested control flow) ----
def fasta(name, text):
    p = os.path.join(tmp, name).replace("\\", "/")
    with open(p, "w", newline="\n") as f:
        f.write(text)
    return p

f1 = fasta("a.fa", ">r1 d\nacgt\nACGT\n\n>r2\nTTTT\n")
f2 = fasta("b.fa", ">x\r\nACGT\r\nNNRY\r\n")
f3 = fasta("c.fa", ">x\nACGJ\n")
f4 = fasta("d.fa", ">x\n")
f5 = fasta("e.fa", ">x\nACGT\nACGT")
f6 = fasta("f.fa", ";comment\n>x\nAC GT\nTT\n")
for nm, p in [("lowercase+2 records", f1), ("crlf+iupac", f2), ("bad char", f3), ("empty", f4),
              ("no trailing newline", f5), ("comment+spaces", f6), ("missing", "/nonexistent/z.fa")]:
    check("fasta " + nm, f'dna g = load_fasta("{p}")\nprint(g)\nprint(len(g))\n')
check("fasta via variable", f'string p = "{f1}"\ndna g = load_fasta(p)\nprint(len(g))\n')

EXTRA = [
    ("shadowing across sibling blocks",
     'int x = 1\nif x == 1 {\nint y = 5\nprint(y)\n}\nif x == 1 {\nfloat y = 2.5\nprint(y)\n}\n'),
    ("loop var scoped", 'dna g = "AAA"\nfor i in find_all(/A/ in g) {\nint z = i * 2\nprint(z)\n}\nfor i in find_all(/A/ in g) {\nprint(i)\n}\n'),
    ("nested loops", 'dna g = "ACACAC"\nfor i in find_all(/A/ in g) {\nfor j in find_all(/C/ in g) {\nif j > i {\nprint(i + j)\n}\n}\n}\n'),
    ("mixed int float", 'int a = 7\nfloat b = 2\nprint(a / b)\nprint(a * 1.5)\nprint(a - 0.5)\nfloat c = a\nprint(c)\n'),
    ("float formats", 'print(1.0 / 3.0)\nprint(1000000.0 * 10.0)\nprint(0.1 + 0.2)\nprint(100.0)\nprint(0.000012345)\n'),
    ("negative via subtraction", 'int a = 3\nprint(0 - a)\nprint(a - 10)\n'),
    ("int division truncation", 'print(7 / 2)\nprint((0 - 7) / 2)\nprint(7 / (0 - 2))\n'),
    ("comparison mixed", 'print(1 < 1.5)\nprint(2 == 2.0)\nprint(3 != 3.0)\nprint(2.5 >= 2)\n'),
    ("bool var", 'bool b = 1 < 2\nprint(b)\nbool c = b == 1 > 2\nprint(c)\n'),
    ("bool compare", 'bool a = 1 < 2\nbool b = 2 < 1\nprint(a == b)\nprint(a != b)\n'),
    ("string ops", 'string a = "foo"\nstring b = "bar"\nprint(a + b)\nprint(len(a + b))\nprint(a == b)\nprint((a + b)[1:4])\n'),
    ("seq eq different length", 'dna a = "AC"\ndna b = "ACG"\nprint(a == b)\n'),
    ("assign reuse", 'dna d = "ACGT"\nd = d + d\nd = d[2:6]\nprint(d)\n'),
    ("protein ops", 'protein p = "MKV*"\nprint(len(p))\nprint(p[1:3])\nprint(find_all(/K/ in p))\n'),
    ("rna find", 'rna r = "AUGAUG"\nprint(find_all(/AUG/ in r))\n'),
    ("motif var print", 'motif m = /A(C|G)+T/\nprint(m)\n'),
    ("empty list loop", 'dna g = "CCC"\nfor i in find_all(/A/ in g) {\nprint(i)\n}\nprint(count(find_all(/A/ in g)))\n'),
    ("empty seq", 'dna g = ""\nprint(find_all(/A/ in g))\nprint(len(g))\nprint(gc_content(g))\nprint(translate(transcribe(g)))\n'),
    ("match_end lazy vs greedy", 'dna g = "ATGAAATAGTAA"\nprint(match_end(/ATG(...)*?(TAA|TAG)/ in g, 0))\nprint(match_end(/ATG(...)*(TAA|TAG)/ in g, 0))\n'),
    ("match_end at end", 'dna g = "ACGT"\nprint(match_end(/T/ in g, 4))\nprint(match_end(/T/ in g, 3))\n'),
    ("extract orf pipeline",
     'dna g = "CCATGAAATTTGGGTAGCC"\nmotif orf = /ATG(<sense>)+<stop>/\nfor i in find_all(orf in g) {\nint e = match_end(orf in g, i)\nprint(i)\nprint(e)\nprint(translate(transcribe(g[i:e])))\n}\n'),
    ("codon classes", 'dna g = "ATGGCCTAAATGTAG"\nprint(find_all(/<codon>{2}<stop>/ in g))\nprint(find_all(/ATG(<sense>)*<stop>/ in g))\n'),
    ("slice in condition", 'dna g = "ACGTACGT"\nif g[0:4] == g[4:8] {\nprint(1)\n} else {\nprint(2)\n}\n'),
    ("else if chain", 'int x = 7\nif x < 3 {\nprint(1)\n} else if x < 5 {\nprint(2)\n} else if x < 9 {\nprint(3)\n} else {\nprint(4)\n}\n'),
    ("err div zero int", 'int a = 0\nprint(5 / a)\n'),
    ("err div zero float", 'float a = 0.0\nprint(5.0 / a)\n'),
    ("err output before error", 'print(1)\nprint(2)\nprint(1 / 0)\nprint(3)\n'),
    ("err slice", 'dna d = "ACGT"\nprint(d[2:9])\n'),
    ("err slice neg", 'dna d = "ACGT"\nint a = 0 - 1\nprint(d[a:2])\n'),
    ("err index", 'dna d = "ACGT"\nprint(find_all(/A/ in d)[3])\n'),
    ("err index neg", 'dna d = "ACGT"\nint a = 0 - 1\nprint(find_all(/A/ in d)[a])\n'),
    ("err match_end pos", 'dna d = "ACGT"\nprint(match_end(/A/ in d, 9))\n'),
    ("err typecheck", 'dna d = "ACGX"\nprint(d)\n'),
    ("err parse", 'int x =\n'),
]
for name, src in EXTRA:
    check("extra: " + name, src)

# ---- 3. motif patterns (same list as motif_check.py) on a 200 kb random sequence ----
PATTERNS = [r"GAG", r"A[CG]T", r"GC+A", r"AT.{3}GC", r"(TAA|TAG|TGA)", r"CG{2,4}T", r"[^AT]{3}", r"TA(T|C)+G",
            r"ATG(...){2,}?(TAA|TAG|TGA)", r"ATG(...){10,}?(TAA|TAG|TGA)", r"ATG(<sense>){2,}<stop>",
            r"ATG(<sense>){10,}<stop>", r"<codon>{4}<stop>", r"A(<sense>)*<stop>", r"[CG]{8}", r"GC.{10}GC"]
random.seed(42)
seq = "".join(random.choice("ACGT") for _ in range(200_000))
for pat in PATTERNS:
    check(f"pattern /{pat}/ on 200 kb",
          f'dna g = "{seq}"\nmotif m = /{pat}/\nprint(find_all(m in g))\nprint(count(find_all(m in g)))\n')

# ---- 4. adversarial inputs: restart budget exceeded -> reverse-DFA linear scan ----
ADV = [("(A+)+B", "A" * 200_000), ("(A|AA)+B", "A" * 100_000), ("A{300}B", "A" * 100_000),
       ("(A*)*C", "A" * 50_000 + "C"), ("(AC|A)*G", "AC" * 40_000 + "A")]
for pat, s in ADV:
    check(f"adversarial /{pat}/ on {len(s)} bases",
          f'dna g = "{s}"\nmotif m = /{pat}/\nprint(find_all(m in g))\n')

# ---- 5. random-pattern fuzzing ----
def rand_pat(depth=0):
    parts = []
    for _ in range(random.randint(1, 4)):
        r = random.random()
        if r < 0.45:
            atom = random.choice("ACGT")
        elif r < 0.6:
            atom = "[" + "".join(random.sample("ACGT", random.randint(1, 3))) + "]"
        elif r < 0.7:
            atom = "."
        elif depth < 2:
            atom = "(" + rand_pat(depth + 1) + ("|" + rand_pat(depth + 1) if random.random() < 0.4 else "") + ")"
        else:
            atom = random.choice("ACGT")
        q = random.random()
        if q < 0.2: atom += "*"
        elif q < 0.35: atom += "+"
        elif q < 0.45: atom += "?"
        elif q < 0.55:
            lo = random.randint(0, 3)
            atom += "{%d,%d}" % (lo, lo + random.randint(0, 3))
        parts.append(atom)
    return "".join(parts)

random.seed(2024)
FUZZ = 150
for k in range(FUZZ):
    pat = rand_pat()
    if random.random() < 0.15:
        pat += "?"                                         # lazy motif
    s = "".join(random.choice("ACGT") for _ in range(random.randint(0, 400)))
    check(f"fuzz /{pat}/ on {len(s)} bases",
          f'dna g = "{s}"\nmotif m = /{pat}/\nprint(find_all(m in g))\nif len(g) > 0 {{\nprint(match_end(m in g, 0))\n}}\n')

print(f"\n{passed} passed, {failed} failed")
print("ALL PASS" if not failed else "SOME FAILED")
sys.exit(1 if failed else 0)
