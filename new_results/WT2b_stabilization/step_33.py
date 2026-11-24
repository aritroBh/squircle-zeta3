import sys, pickle
sys.path.insert(0,'.')
from fastred import *
record = pickle.load(open('record_W2b.pkl','rb'))

def seq_H(key, r):
    d=sum(key); out={}
    for k in range(d, 13):
        rec = record.get((k,key))
        if rec is None: continue
        H = rec['H'][1]
        out[k] = H[r-1] if len(H)>=r else ZERO
    return out

def min_poly(seq, maxdeg):
    # minimal (k0, deg<=maxdeg) s.t. tail from k0 is poly deg<=deg, verified with >=deg+2 pts
    ks = sorted(seq)
    best=None
    for k0 in ks:
        kk=[k for k in ks if k>=k0]; vals=[seq[k] for k in kk]
        for deg in range(-1, maxdeg+1):
            need = deg+2 if deg>=0 else 1
            if len(vals) < max(need,2): continue
            ds = vals
            for _ in range(max(deg,0)):
                ds=[ds[t+1]-ds[t] for t in range(len(ds)-1)]
            if deg==-1:
                ok = all(v.is_zero() for v in ds)
            else:
                ok = all((ds[t+1]-ds[t]).is_zero() for t in range(len(ds)-1))
            if ok:
                if best is None or k0<best[0] or (k0==best[0] and deg<best[1]):
                    best=(k0,deg)
                break
    return best

print("H (higher part of W_ij): polynomial behavior of +i coefficients")
for key in [(2,0),(1,1),(3,0),(2,1),(4,0),(3,1),(2,2),(5,0),(4,1),(3,2),(6,0),(5,1),(4,2),(3,3)]:
    d=sum(key); line=[]
    for r in range(2, 2*d-1):
        bound = d-1-(-(-r)//2)
        b = min_poly(seq_H(key,r), bound+1)
        line.append(f"r{r}:k0={b[0]},deg{b[1]}/B{bound}" + ('' if b[1]<=bound and b[0]<=d+1 else '!'))
    print(f"{key}: "+"  ".join(line))
