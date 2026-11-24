
import sys, os, json, glob, time
print("python:", sys.executable)
print(sys.version)
import mpmath as mp
print("mpmath", mp.__version__)
try:
    import flint
    print("flint", flint.__version__)
    print("lll sig/doc:", (flint.fmpz_mat.lll.__doc__ or "")[:400])
except Exception as e:
    print("flint import error", repr(e))
import sympy
print("sympy", sympy.__version__)
print("cwd:", os.getcwd())
print("top listing:", sorted(os.listdir("."))[:60])
# look for the referenced files, depth-limited
targets = []
for base in ["."]:
    for root, dirs, files in os.walk(base):
        depth = root.count(os.sep)
        if depth > 6:
            dirs[:] = []
            continue
        for f in files:
            if f in ("L_relation.json",) or (f=="result.md" and "M4i" in root) or f=="report_b.md":
                targets.append(os.path.join(root, f))
print("targets found:")
for t in targets[:40]:
    print("  ", t)
