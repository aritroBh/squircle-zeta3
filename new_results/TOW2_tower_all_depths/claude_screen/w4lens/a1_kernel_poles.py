# Item (2): +-i simple pole / const content of degree-3 (and 4,5) kernels, all available k.
import pickle, os, sys
sys.dont_write_bytecode = True
TH="<project>/research/08_theory/WTS_weight_k_structure/theory"
RD="<project>/research/08_theory/TWOEND1_weight_poles/red"
def path(k):
    p=os.path.join(TH,f"kern_{k}.pkl")
    if os.path.exists(p): return p
    p=os.path.join(RD,f"kern_{k}.pkl")
    return p if os.path.exists(p) else None
def nz(t): return any(x!=0 for x in t)
avail=[k for k in range(1,30) if path(k)]
print("levels available:",avail)
for d in (1,2,3,4,5):
  print("== degree",d)
  for k in avail:
    K=pickle.load(open(path(k),'rb'))
    keys=sorted([key for key in K if sum(key)==d])
    row=[]
    for key in keys:
        e=K[key]
        ip=[nz(t) for t in e['ipos']]; ig=[nz(t) for t in e['ineg']]
        row.append((key, 'ipos_nonzero' if any(ip) else 'ipos0', 'ineg_nonzero' if any(ig) else 'ineg0', 'const_nz' if nz(e['const']) else 'const0', 'len(ipos)=%d'%len(e['ipos'])))
    print(f"k={k:2d}", row)
