
import os, pickle
BASE = '<project>'
WT2B = os.path.join(BASE, 'research/08_theory/WT2b_stabilization')
rec = pickle.load(open(os.path.join(WT2B,'record_W2b.pkl'),'rb'))

for key in [(3,0),(2,1)]:
    print(f"===== Gm for {key} =====")
    ks = sorted(k for (k,kk) in rec if kk==key)
    rows = {}
    for k in ks:
        G = rec[(k,key)]['Gm']
        for m, c in enumerate(G[1], start=1):
            rows.setdefault(('ipos',m), {})[k] = c
        rows.setdefault(('const',0), {})[k] = G[3]
    for tag in sorted(rows):
        print(" ", tag, {k: str(c) for k,c in sorted(rows[tag].items())})
