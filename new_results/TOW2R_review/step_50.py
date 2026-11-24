import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
os.chdir(os.path.join(workdir, 'tow2_copy'))
from fractions import Fraction as Fr
import genintegrand as GI

# Fr-polynomial arithmetic (low-to-high)
def padd(P,Q):
    n=max(len(P),len(Q)); out=[Fr(0)]*n
    for k,c in enumerate(P): out[k]+=c
    for k,c in enumerate(Q): out[k]+=c
    return out
def pmul(P,Q):
    if not P or not Q: return []
    out=[Fr(0)]*(len(P)+len(Q)-1)
    for a,ca in enumerate(P):
        if ca:
            for b,cb in enumerate(Q):
                if cb: out[a+b]+=ca*cb
    return out
def pscale(P,q): return [c*q for c in P]
def psub(P,Q): return padd(P, pscale(Q, Fr(-1)))

# v-polynomials
V = [Fr(0), Fr(1)]
def vpow(P, n):
    out=[Fr(1)]
    for _ in range(n): out=pmul(out,P)
    return out

# t(v) = (v^2+2v-1)/(4v);  Stilde(v) = sum_j s_j (v^2+2v-1)^j (4v)^{D-j}, D=deg S
def Stilde(S):
    D = len(S)-1
    base = [Fr(-1), Fr(2), Fr(1)]     # v^2+2v-1
    fourv = [Fr(0), Fr(4)]
    out = []
    for j, c in enumerate(S):
        if c==0: continue
        term = pmul(vpow(base, j), vpows(fourv, D-j)) if False else pmul(vpow(base,j), [Fr(4)**(D-j)] + [Fr(0)]*(0))
        # (4v)^{D-j} = 4^{D-j} v^{D-j}: shift
        t2 = pscale(vpow(base, j), c * Fr(4)**(D-j))
        t2 = [Fr(0)]*(D-j) + t2
        out = padd(out, t2)
    return out

# check deg S = 2d for all keys, d=2..12
degs = set()
for dd in range(2, 13):
    for mm in range(0, dd+1):
        S = GI.S_poly(dd, mm)
        degs.add(len(S)-1)
print('deg S over d=2..12 all m:', sorted(degs))
