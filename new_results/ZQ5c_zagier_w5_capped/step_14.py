base = '<project>/research/08_theory'
import json, os
# Check d29 and d13 candidates
for dd in ['d29','d13','d21']:
    p = f'{base}/ZQ5b_zagier_w5_larger_S/out/{dd}/candidates.json'
    J = json.load(open(p))
    print("="*60)
    print(dd, "keys:", list(J.keys()))
    print("names:", J.get('names'))
    print("n cands:", len(J.get('cands',[])))
    for c in J.get('cands',[])[:5]:
        print("  ", c)
# check d17_final format
print("\n===== d17_final.json (head) =====")
J = json.load(open(f'{base}/ZQ5b_zagier_w5_larger_S/d17_final.json'))
print("keys:", list(J.keys()))
for k in J:
    v = J[k]
    if isinstance(v, dict):
        print(f"  {k}: dict with {len(v)} entries; first 3:", dict(list(v.items())[:3]))
    else:
        print(f"  {k}: {v}")
