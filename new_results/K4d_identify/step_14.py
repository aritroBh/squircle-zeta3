import re
SRC = 'research/08_theory/K4c_3d_edge_law/'
tr = open(SRC+'transcript.json').read()
# model B d^k values: search for leading digits
for pat in ['3.5164954', '2.3544806', '1.8983271', '0.4126500']:
    idxs=[m.start() for m in re.finditer(re.escape(pat), tr)]
    print(f"{pat}: {len(idxs)} occ")
    for i in idxs[:2]:
        print("   ", tr[i:i+75])
