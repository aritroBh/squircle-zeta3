
import os, json
base = "<project>"
h1 = json.load(open(os.path.join(base,"research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json")))
print(h1.keys())
for k,v in h1.items():
    if isinstance(v, list):
        print(k, "len", len(v), "first:", str(v[0])[:60], "last:", str(v[-1])[:80])
    else:
        print(k, "=", str(v)[:200])
