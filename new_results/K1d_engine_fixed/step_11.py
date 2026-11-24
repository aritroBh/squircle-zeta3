import json
for jname in ['t3_b0_150.json','t3_b1_110.json']:
    d = json.load(open(jname))
    print("="*20, jname)
    print("top-level keys:", list(d.keys()))
    arc = d['arcs'][0]
    print("arc keys:", list(arc.keys()))
    for k,v in arc.items():
        s = str(v)
        print(f"  {k}: {s[:120]}")
