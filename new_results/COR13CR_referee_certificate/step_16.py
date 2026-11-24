import os, json
base = '<project>/research/08_theory'
d = os.path.join(base,'COR13C_third_corner_certificate')
c23 = json.load(open(os.path.join(d,'cert13_w23.json')))
print("cert23 keys:", list(c23.keys()))
for k in c23:
    C = c23[k]
    print(k, "keys:", list(C.keys()), "prow len:", len(C['prow']), "y len:", len(C['y']), "cols len:", len(C['cols']), "target support:", len(C.get('target',[])))
c4 = json.load(open(os.path.join(d,'cert13_w4.json')))
print("cert4 keys:", list(c4.keys()))
print("prow:", len(c4['prow']), "y:", len(c4['y']), "cols:", len(c4['cols']), "target:", len(c4.get('target',[])))
print("y sample:", c4['y'][:3])
euler = json.load(open(os.path.join(d,'cert13_euler.json')))
print("euler:", {k: (len(v) if isinstance(v,list) else v) for k,v in euler.items()})
mc = json.load(open(os.path.join(d,'map_counts.json')))
print("map_counts:", str(mc)[:600])
