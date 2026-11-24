# step_30: inspect 12-gon CR set, orbits of the 8 DOD5 points, duplication relations
from fractions import Fraction as Fr
import pickle
CR12 = pickle.load(open("CR_regular.pkl","rb"))["CR12"]
import sys
sys.path.insert(0,".")
# rebuild QF F3 quickly
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp
def qmul(x, y, D): return (x[0]*y[0] + D*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def qpow(x, n, D):
    if n < 0:
        a,b = x; den = a*a - D*b*b
        return qpow((a/den, -b/den), -n, D)
    r = (Fr(1), Fr(0))
    while n:
        if n & 1: r = qmul(r, x, D)
        x = qmul(x, x, D); n >>= 1
    return r
def mpf_of(q): q=Fr(q); return mp.mpf(q.numerator)/mp.mpf(q.denominator)
def vp(q, p):
    q = abs(Fr(q)); n, d = q.numerator, q.denominator
    e = 0
    while n % p == 0 and n>1: n//=p; e+=1
    while d % p == 0 and d>1: d//=p; e-=1
    return e
class QF:
    def __init__(self, D, letters):
        self.D = Fr(D); self.letters = letters
        self.LAB = {n: ab for n, ab, _ in letters}
        self.LNM = {n: nrm for n, _, nrm in letters}
    def elem(self, v):
        r = (Fr(1), Fr(0))
        for nm, e in v.items(): r = qmul(r, qpow(self.LAB[nm], int(e), self.D), self.D)
        return r
    def vec(self, z):
        a,b = Fr(z[0]), Fr(z[1])
        N = a*a - self.D*b*b
        exps = {}
        for nm, ab, nrm in self.letters:
            if abs(nrm) == 1: continue
            fct = sp.factorint(abs(nrm)); pr, k = next(iter(fct.items()))
            v = vp(N, pr); assert v % k == 0
            exps[nm] = v//k
        rem = abs(N)
        for nm, e in exps.items(): rem /= abs(self.LNM[nm])**e
        assert rem == 1, (z, rem)
        mp.mp.dps = 40
        sq = mp.sqrt(int(self.D))
        def aval(ab): return abs(mpf_of(ab[0]) + sq*mpf_of(ab[1]))
        lt = mp.log(aval((a,b)))
        prod_log = mp.mpf(0)
        for nm, ab, nrm in self.letters:
            e = exps.get(nm, 0)
            if e and abs(nrm) != 1: prod_log += e*mp.log(aval(ab))
        units = [(nm,ab) for nm,ab,nrm in self.letters if abs(nrm)==1]
        resid = lt - prod_log
        if units:
            nm, ab = units[0]
            c = resid/mp.log(aval(ab)); ci = int(mp.nint(c))
            assert abs(c-ci) < mp.mpf('1e-25'), (z,c)
            exps[nm] = exps.get(nm,0)+ci
        got = self.elem(exps)
        assert got == (a,b) or got == (-a,-b), (z, exps, got)
        return {k:Fr(v) for k,v in exps.items() if v}
F3 = QF(3, [("g2",(Fr(1),Fr(1)),-2), ("g3",(Fr(0),Fr(1)),-3), ("u",(Fr(2),Fr(1)),1)])

def qsub(x,y,D): return (x[0]-y[0], x[1]-y[1])
def qinv(x,D):
    a,b=x; den=a*a-D*b*b; return (a/den,-b/den)
def qdiv(x,y,D): return qmul(x,qinv(y,D),D)
def S3orbit(z, D):
    pts = [z, qsub((Fr(1),Fr(0)),z,D)]
    if z != (0,0): pts.append(qinv(z,D))
    if z != (1,0): pts.append(qinv(qsub((Fr(1),Fr(0)),z,D),D))
    if z not in [(0,0),(1,0)]:
        pts += [qdiv(z,qsub(z,(Fr(1),Fr(0)),D),D), qdiv(qsub(z,(Fr(1),Fr(0)),D),z,D)]
    return set(pts)

print("CR_12 with vecs:")
for z in sorted(CR12):
    print("   ", z, F3.vec(z))
z8 = [(2,-1),(Fr(-1,2),Fr(1,2)),(0,Fr(1,3)),(7,-4),(1,Fr(-1,2)),(-1,Fr(2,3)),(Fr(1,2),Fr(-1,6)),(Fr(1,3),0)]
z8 = [tuple(Fr(t) for t in z) for z in z8]
names = ["2-s3","(s3-1)/2","s3/3","7-4s3","(2-s3)/2","(2s3-3)/3","(3-s3)/6","1/3"]
print("\nDOD5 points: orbit intersection with CR, duplication hits:")
for nm,z in zip(names,z8):
    orb = S3orbit(z, 3)
    inter = orb & CR12
    # duplication: z = w^2 or -w^2 or z^2 in CR
    dups = []
    for w in CR12:
        w2 = qmul(w,w,3)
        if w2 == z: dups.append(('sqrt', w))
        if qmul(w2,(-1,0),3) == z or (w2[0]==-z[0] and w2[1]==-z[1]): dups.append(('neg sqrt', w))
    z2 = qmul(z,z,3)
    if z2 in CR12: dups.append(('square in CR', z2))
    print(f"  {nm} {z}: orbit cap CR = {sorted(inter)}; dup: {dups}")
