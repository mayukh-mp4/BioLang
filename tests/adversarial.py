"""Adversarial inputs: Python's backtracking re versus the BioLang DFA engine.

Usage (from the project root):
    python tests\\adversarial.py build\\biolang.exe

Case 1 is the classic catastrophic-backtracking pattern (A+)+B on a string with no B.
Case 2 is a long partial match that defeats naive restart scanning: A{N}B on a long run of A's.
Times are matching time only, in milliseconds.
"""
import os, re, subprocess, sys, tempfile, time

exe = sys.argv[1] if len(sys.argv) > 1 else ("build/biolang.exe" if os.name == "nt" else "build/biolang")
PY_TIMEOUT_N = 24          # length of the run of A's given to Python in case 1 (it is exponential)
tmp = tempfile.mkdtemp()


def run_biolang(seq, pat, mode):
    path = os.path.join(tmp, "adv.bl")
    with open(path, "w", newline="\n") as f:
        f.write(f'dna s = "{seq}"\nmotif m = /{pat}/\nprint(count(find_all(m in s)))\n')
    r = subprocess.run([exe, "--time", f"--scan={mode}", path], capture_output=True, text=True)
    line = [l for l in r.stderr.splitlines() if l.startswith("[time]")][0]
    t = dict(kv.split("=") for kv in line.replace("[time]", "").split())
    return float(t["scan_ms"]), int(r.stdout.strip())


def run_python(seq, pat):
    rx = re.compile(pat)
    t = time.perf_counter()
    n = sum(1 for _ in rx.finditer(seq))
    return (time.perf_counter() - t) * 1000, n


def row(label, py, restart, auto, linear):
    f = lambda v: "n/a" if v is None else f"{v:,.1f} ms"
    print(f"{label:<34} {f(py):>14} {f(restart):>14} {f(auto):>12} {f(linear):>12}")


print(f"{'case':<34} {'Python re':>14} {'BioLang restart':>14} {'BioLang auto':>12} {'BioLang linear':>12}")

# Case 1: (A+)+B  on A...AC  -> exponential for backtracking engines
pat = "(A+)+B"
small = "A" * PY_TIMEOUT_N + "C"
py, n1 = run_python(small, pat)
rs, n2 = run_biolang(small, pat, "restart"); au, _ = run_biolang(small, pat, "auto"); li, _ = run_biolang(small, pat, "linear")
assert n1 == n2 == 0
row(f"(A+)+B, {PY_TIMEOUT_N} A's then C", py, rs, au, li)
big = "A" * 200_000 + "C"
rs, _ = run_biolang(big, pat, "restart"); au, _ = run_biolang(big, pat, "auto"); li, _ = run_biolang(big, pat, "linear")
row("(A+)+B, 200,000 A's then C", None, rs, au, li)

# Case 2: A{N}B on a long run of A's -> every start position has a long partial match
for N in (300, 1000):
    pat = "A{%d}B" % N
    seq = "A" * 200_000
    py, n1 = run_python(seq, pat)
    rs, n2 = run_biolang(seq, pat, "restart"); au, n3 = run_biolang(seq, pat, "auto"); li, n4 = run_biolang(seq, pat, "linear")
    assert n1 == n2 == n3 == n4 == 0
    row(f"A{{{N}}}B on 200,000 A's", py, rs, au, li)
