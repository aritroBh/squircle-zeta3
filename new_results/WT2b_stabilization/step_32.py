import sys, pickle
sys.path.insert(0,'.')
from fastred import *
record = pickle.load(open('record_W2b.pkl','rb'))

# collect sequences: for each canonical key (i,j) and pole order r>=2 of H (higher part)
# sequence over k >= i+j of the +i coefficient
from fractions import Fraction as Fr

def seq_H(key, r):
    d=sum(key); out={}
    for k in range(d, 13):
        rec = record.get((k,key))
        if rec is None: continue
        H = rec['H'][1]
        out[k] = H[r-1] if len(H)>=r else ZERO
    return out

def seq_G(key, s):
    d=sum(key); out={}
    for k in range(d, 13):
        rec = record.get((k,key))
        if rec is None or 'Gm' not in rec: continue
        G = rec['Gm'][1]
        out[k] = G[s-1] if len(G)>=s else ZERO
    return out

def fit_poly(seq, maxdeg):
    # find min deg <= maxdeg and min start k0 such that seq[k] for k>=k0 is poly of that deg
    ks = sorted(seq)
    for deg in range(-1, maxdeg+1):
        # try start offsets
        for k0 in ks:
            vals = [seq[k] for k in ks if k>=k0]
            kk = [k for k in ks if k>=k0]
            if deg==-1:
                if all(v.is_zero() for v in vals): return (-1,k0)
                continue
            if len(vals) < deg+2:  # need deg+2 points to confirm
                # still try interpolation+check on available
                pass
            # finite differences: deg-th diff constant
            ds = vals
            ok=True
            for _ in range(deg):
                ds = [ds[t+1]-ds[t] for t in range(len(ds)-1)]
            if ds and all((ds[t+1]-ds[t]).is_zero() for t in range(len(ds)-1)) and len(ds)>=1:
                return (deg,k0)
    return (None,None)

print("KEY  ord  bound  (fitdeg, k0)   [H = higher part of W]")
for key in [(2,0),(1,1),(3,0),(2,1),(4,0),(3,1),(2,2),(5,0),(4,1),(3,2),(6,0),(5,1),(4,2),(3,3)]:
    d=sum(key)
    for r in range(2, 2*d-1):
        bound = d-1-(-(-r)//2)
        sq = seq_H(key,r)
        if not sq: continue
        fd,k0 = fit_poly(sq, bound+1)
        flag = '' if (fd is not None and fd<=bound and k0<=d+1) else '  <<< CHECK'
        print(f"{key} r={r} B={bound}: fit deg {fd} from k0={k0}{flag}")
