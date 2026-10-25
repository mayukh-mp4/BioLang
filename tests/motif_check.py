"""Correctness check: BioLang's DFA motif engine vs Python's `re`.

Usage (from the project root):
    python tests\\motif_check.py build\\biolang.exe
"""
import os, random, re, subprocess, sys, tempfile

exe = sys.argv[1] if len(sys.argv) > 1 else ("build/biolang.exe" if os.name == "nt" else "build/biolang")

# Patterns chosen so that Python's backtracking rules and BioLang's DFA rules
# (leftmost, longest - or shortest for lazy) give the same match.
PATTERNS = [
    r"GAG",
    r"A[CG]T",
    r"GC+A",
    r"AT.{3}GC",
    r"(TAA|TAG|TGA)",
    r"CG{2,4}T",
    r"[^AT]{3}",
    r"TA(T|C)+G",
    r"ATG(...){2,}?(TAA|TAG|TGA)",
    r"ATG(...){10,}?(TAA|TAG|TGA)",
]

random.seed(42)
seq = "".join(random.choice("ACGT") for _ in range(200_000))
tmp = tempfile.mkdtemp()
failed = 0
for pat in PATTERNS:
    path = os.path.join(tmp, "t.bl")
    with open(path, "w", newline="\n") as f:
        f.write(f'dna g = "{seq}"\nmotif m = /{pat}/\nprint(find_all(m in g))\n')
    res = subprocess.run([exe, path], capture_output=True, text=True)
    got = res.stdout.strip()
    want = "[" + ", ".join(str(m.start()) for m in re.finditer(pat, seq) if m.end() > m.start()) + "]"
    ok = got == want
    failed += not ok
    n = want.count(",") + 1 if want != "[]" else 0
    print(f"{'PASS' if ok else 'FAIL'}  /{pat}/  ({n} matches)")
    if not ok:
        print("   stderr:", res.stderr.strip()[:200])
print("\nAll patterns agree with Python re." if not failed else f"\n{failed} pattern(s) differ.")
sys.exit(1 if failed else 0)
