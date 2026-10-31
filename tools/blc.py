"""BioLang -> native executable via LLVM.

Usage (from the project root):
    python tools\\blc.py program.bl [-o program.exe] [--keep-ll] [--run]

Needs clang (LLVM >= 15) on PATH; set CLANG to use another binary.
Steps: biolang --emit-llvm program.bl > program.ll ; clang -O2 program.ll runtime/biolang_rt.c
"""
import argparse, os, subprocess, sys

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
default_exe = os.path.join(root, "build", "biolang.exe" if os.name == "nt" else "biolang")

ap = argparse.ArgumentParser()
ap.add_argument("source")
ap.add_argument("-o", dest="out")
ap.add_argument("--biolang", default=default_exe)
ap.add_argument("--keep-ll", action="store_true")
ap.add_argument("--run", action="store_true")
a = ap.parse_args()

base = os.path.splitext(a.source)[0]
out = a.out or base + (".exe" if os.name == "nt" else ".out")
ll = base + ".ll"
clang = os.environ.get("CLANG", "clang")

r = subprocess.run([a.biolang, "--emit-llvm", a.source], capture_output=True, text=True)
if r.returncode != 0:
    sys.stderr.write(r.stderr)
    sys.exit(1)
with open(ll, "w", newline="\n") as f:
    f.write(r.stdout)
c = subprocess.run([clang, "-O2", "-Wno-override-module", ll, os.path.join(root, "runtime", "biolang_rt.c"), "-o", out])
if not a.keep_ll:
    os.remove(ll)
if c.returncode != 0:
    sys.exit(c.returncode)
print(f"wrote {out}")
if a.run:
    sys.exit(subprocess.run([out]).returncode)
