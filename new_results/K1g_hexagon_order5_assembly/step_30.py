import os, json
base = '<project>/research/08_theory'
k1d = os.path.join(base, 'K1d_engine_fixed')
print("K1d files:", sorted(os.listdir(k1d)))
print()
# look at geometry.json
p = os.path.join(k1d, 'geometry.json')
if os.path.exists(p):
    with open(p) as f:
        g = json.load(f)
    print("geometry.json keys:", list(g.keys()))
    import pprint
    pprint.pprint(g[list(g.keys())[0]][0])
