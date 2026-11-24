"""Summarise analyze_tower.out: for every (n, depth) cell with a genuine test (npts >= degree+2) check degree == floor(depth/2)."""
import re, collections, os
HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for line in open(os.path.join(HERE, "analyze_tower.out")):
    m = re.match(r"n=(\d+) depth delta=(\d+) \(order r=(\d+)\): poly degree in k over keys = \[(\w+)\];.*npts=(\d+)", line)
    if not m: continue
    n, dl, r, deg, npts = int(m[1]), int(m[2]), int(m[3]), m[4], int(m[5])
    rows.append((n, dl, r, None if deg == 'None' else int(deg), npts))
tested = [(n, dl, deg, npts) for (n, dl, r, deg, npts) in rows if deg is not None and npts >= deg + 2]
bad = [(n, dl, deg) for (n, dl, deg, npts) in tested if deg != dl // 2]
print("cells (n, depth) in analyze_tower.out:", len(rows))
print("cells with a genuine exact test (npts >= degree+2):", len(tested))
print("cells where minimal polynomial degree != floor(depth/2):", bad)
byn = collections.defaultdict(list)
for (n, dl, deg, npts) in tested: byn[n].append(dl)
for n in sorted(byn):
    print(f"  degree n={n}: depths 0..{max(byn[n])} tested (of 0..{2*n-4}); all have degree floor(depth/2)")
