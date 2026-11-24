
import os, pickle
BASE = '<project>'
WT2B = os.path.join(BASE, 'research/08_theory/WT2b_stabilization')
rec = pickle.load(open(os.path.join(WT2B,'record_W2b.pkl'),'rb'))

def show(key, what='H'):
    print(f"==== {what} for key {key}, +i principal parts by order ====")
    # collect per order across k
    per_order = {}
    ks = sorted(k for (k, kk) in rec if kk == key)
    for k in ks:
        P = rec[(k, key)][what]
        ipos = P[1]
        for m, c in enumerate(ipos, start=1):
            per_order.setdefault(m, {})[k] = c
    for m in sorted(per_order):
        row = per_order[m]
        print(f" order {m}: ", {k: str(c) for k, c in sorted(row.items())})

show((3,0), 'H')
show((2,1), 'H')
