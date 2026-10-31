"""Fuzz test: the linear scan must give exactly the same matches as the restart scan.

Usage (from the project root):
    python tests\\scan_equivalence.py build\\biolang.exe

Random motifs (with alternation, groups, classes, counted and lazy quantifiers) are run over random
and low-complexity DNA, once with --scan=linear and once with --scan=restart; the outputs must match.
"""
import os, random, subprocess, sys, tempfile

exe = sys.argv[1] if len(sys.argv) > 1 else ("build/biolang.exe" if os.name == "nt" else "build/biolang")
random.seed(7)


def atom(depth):
    r = random.random()
    if r < 0.45: return random.choice("ACGT")
    if r < 0.55: return "."
    if r < 0.70: return "[" + "".join(random.sample("ACGT", random.randint(1, 3))) + "]"
    if r < 0.75: return "[^" + random.choice("ACGT") + "]"
    if depth < 2: return "(" + alt(depth + 1) + ")"
    return random.choice("ACGT")


def quant(a):
    if random.random() < 0.5: return a
    q = random.choice(["*", "+", "?", "{2}", "{1,3}", "{2,}"])
    return a + q + ("?" if random.random() < 0.15 else "")


def cat(depth): return "".join(quant(atom(depth)) for _ in range(random.randint(1, 4)))
def alt(depth): return "|".join(cat(depth) for _ in range(random.choice([1, 1, 2, 3])))


def random_seq():
    kind = random.random()
    if kind < 0.5: return "".join(random.choice("ACGT") for _ in range(300))
    if kind < 0.8: return "".join(random.choice("AAAC") for _ in range(300))   # low complexity
    unit = "".join(random.choice("ACGT") for _ in range(random.randint(2, 6)))
    return (unit * 100)[:300]


tmp = tempfile.mkdtemp()
cases, bad, skipped = 0, 0, 0
batch = 0
while batch < 40:
    lines = []
    for k in range(3):
        lines.append(f'dna s{k} = "{random_seq()}"')
    for m in range(25):
        pat = alt(0)
        lines.append(f"motif m{m} = /{pat}/")
        for k in range(3):
            lines.append(f"print(find_all(m{m} in s{k}))")
    path = os.path.join(tmp, "fuzz.bl")
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    out, too_complex = {}, False
    for mode in ("linear", "restart"):
        r = subprocess.run([exe, f"--scan={mode}", path], capture_output=True, text=True)
        if r.returncode != 0:
            if "too complex" in r.stderr:      # a random pattern exceeded the DFA size cap: skip the batch
                too_complex = True
                break
            print("ERROR in batch", batch, mode, r.stderr[:200]); sys.exit(2)
        out[mode] = r.stdout.splitlines()
    if too_complex:
        skipped += 1
        continue
    batch += 1
    cases += len(out["restart"])
    for i, (a, b) in enumerate(zip(out["linear"], out["restart"])):
        if a != b:
            bad += 1
            print(f"MISMATCH batch {batch} output line {i}:\n  linear : {a[:100]}\n  restart: {b[:100]}")
print(f"{cases} motif/sequence pairs compared, {bad} mismatches ({skipped} batches skipped: pattern too complex)")
sys.exit(1 if bad else 0)
