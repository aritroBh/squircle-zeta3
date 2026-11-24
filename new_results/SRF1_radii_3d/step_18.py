
import json
h1 = '<project>/research/10_frontier/H1_octahedron_end'
j3 = json.load(open(h1+'/hi_coeffs_n3_K80.json'))
print("keys:", list(j3.keys()))
for k in list(j3.keys()):
    v = j3[k]
    if isinstance(v, list): print(k, "list len", len(v), str(v[:3])[:200])
    else: print(k, str(v)[:120])
