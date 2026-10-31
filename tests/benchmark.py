"""Speed + correctness benchmark: BioLang DFA engine vs Python `re`.

Usage (from the project root):
    python tests\\benchmark.py build\\biolang.exe genome1.fa [genome2.fa ...]
    python tests\\benchmark.py build\\biolang.exe --synthetic 1 5 20      (random genomes, sizes in Mbp)

Both sides time ONLY the matching step:
  * Python: compiled_pattern.finditer(seq) over the whole sequence (re.compile excluded)
  * BioLang: find_all scan time reported by `biolang --time` (file loading and DFA build excluded;
    the DFA build time is shown separately)
Each measurement is the best of REPS runs. Match counts must agree, or the row is flagged.
Results are also written to bench/results.csv.
"""
import csv, os, random, re, subprocess, sys, tempfile, time

REPS = 3

# BioLang codon classes and their plain-regex equivalents (for Python's re)
MACROS = {
    "<codon>": "(?:[ACGT][ACGT][ACGT])",
    "<stop>": "(?:TAA|TAG|TGA)",
    "<sense>": "(?:[ACG][ACGT][ACGT]|T(?:[CT][ACGT]|A[CT]|G[CGT]))",
}


def to_python(pat):
    for k, v in MACROS.items():
        pat = pat.replace(k, v)
    return pat


PATTERNS = [
    r"GAATTC",                                # EcoRI restriction site (plain literal)
    r"TATA[AT]A[AT]",                         # TATA-box-like motif
    r"[CG]{8}",                               # GC-rich run
    r"GC.{10}GC",                             # spaced motif
    r"ATG(...){10,}?(TAA|TAG|TGA)",           # ORF, >= 10 codons
    r"ATG(...){30,}?(TAA|TAG|TGA)",           # ORF, >= 30 codons
    r"ATG(<sense>){30,}<stop>",               # strict ORF: no internal stop codon
]


def read_fasta(path):
    seq = []
    started = False
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if started:
                    break
                started = True
            elif line:
                started = True
                seq.append(line.upper())
    return "".join(seq)


def make_synthetic(mbp, outdir):
    random.seed(mbp)
    seq = "".join(random.choices("ACGT", k=int(mbp * 1_000_000)))
    path = os.path.join(outdir, f"synthetic_{mbp}Mbp.fa")
    with open(path, "w", newline="\n") as f:
        f.write(f">synthetic random genome {mbp} Mbp\n")
        for i in range(0, len(seq), 70):
            f.write(seq[i:i + 70] + "\n")
    return path


def bench_python(pat, seq):
    rx = re.compile(to_python(pat))
    best, count = None, 0
    for _ in range(REPS):
        t = time.perf_counter()
        count = sum(1 for m in rx.finditer(seq) if m.end() > m.start())
        dt = (time.perf_counter() - t) * 1000
        best = dt if best is None else min(best, dt)
    return best, count


def bench_biolang(exe, fasta, pat, tmpdir):
    prog = os.path.join(tmpdir, "bench.bl")
    fa = os.path.abspath(fasta).replace("\\", "/")
    with open(prog, "w", newline="\n") as f:
        f.write(f'dna g = load_fasta("{fa}")\nmotif m = /{pat}/\nprint(count(find_all(m in g)))\n')
    best_scan, best_compile, count = None, None, 0
    for _ in range(REPS):
        res = subprocess.run([exe, "--time", prog], capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(res.stderr.strip())
        count = int(res.stdout.strip())
        line = [l for l in res.stderr.splitlines() if l.startswith("[time]")][0]
        t = dict(kv.split("=") for kv in line.replace("[time]", "").split())
        scan, comp = float(t["scan_ms"]), float(t["compile_ms"])
        best_scan = scan if best_scan is None else min(best_scan, scan)
        best_compile = comp if best_compile is None else min(best_compile, comp)
    return best_scan, best_compile, count


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print(__doc__)
        sys.exit(1)
    exe, rest = args[0], args[1:]
    tmpdir = tempfile.mkdtemp()
    if rest[0] == "--synthetic":
        genomes = [make_synthetic(float(x), tmpdir) for x in rest[1:]]
    else:
        genomes = rest

    rows = []
    for g in genomes:
        seq = read_fasta(g)
        print(f"\nGenome: {os.path.basename(g)}  ({len(seq):,} bp)")
        print(f"{'pattern':34} {'matches':>9} {'Python re':>11} {'BioLang':>10} {'speedup':>8} {'DFA build':>10}")
        for pat in PATTERNS:
            py_ms, py_n = bench_python(pat, seq)
            bl_ms, comp_ms, bl_n = bench_biolang(exe, g, pat, tmpdir)
            flag = "" if py_n == bl_n else f"   MISMATCH (python={py_n}, biolang={bl_n})"
            print(f"/{pat}/".ljust(34) + f" {bl_n:>9,} {py_ms:>9.1f}ms {bl_ms:>8.1f}ms {py_ms / bl_ms:>7.1f}x {comp_ms:>8.2f}ms{flag}")
            rows.append([os.path.basename(g), len(seq), pat, bl_n, round(py_ms, 3), round(bl_ms, 3), round(comp_ms, 3)])

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bench")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "results.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["genome", "bp", "pattern", "matches", "python_ms", "biolang_ms", "dfa_build_ms"])
        w.writerows(rows)
    print(f"\nResults written to {os.path.normpath(out)}")


if __name__ == "__main__":
    main()
