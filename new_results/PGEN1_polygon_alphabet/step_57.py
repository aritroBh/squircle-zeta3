# step_32: 12-gon symbol span check over letters {g2,g3,u}
from fractions import Fraction as Fr
import pickle, itertools
import sympy as sp

CR12 = pickle.load(open("CR_regular.pkl","rb"))["CR12"]
# QF vecfun (exact, verified) -- reuse from step_28
import mpmath as mp
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
        self.D = Fr(D); self.LAB = {n: ab for n, ab, _ in letters}
        self.LNM = {n: nrm for n, _, nrm in letters}; self.letters = letters
    def elem(self, v):
        r = (Fr(1), Fr(0))
        for nm, e in v.items(): r = qmul(r, qpow(self.LAB[nm], int(e), self.D), self.D)
        return r
    def vec(self, z):
        a,b = Fr(z[0]), Fr(z[1]); N = a*a - self.D*b*b
        exps = {}
        for nm, ab, nrm in self.letters:
            if abs(nrm) == 1: continue
            fct = sp.factorint(abs(nrm)); pr, k = next(iter(fct.items()))
            v = vp(N, pr); assert v % k == 0, (z,nm,v)
            exps[nm] = v//k
        rem = abs(N)
        for nm, e in exps.items(): rem /= abs(self.LNM[nm])**e
        assert rem == 1, (z, rem)
        mp.mp.dps = 40
        sq = mp.sqrt(int(self.D))
        def aval(ab): return abs(mpf_of(ab[0]) + sq*mpf_of(ab[1]))
        lt = mp.log(aval((a,b))); prod_log = mp.mpf(0)
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

def qsub(x,y): return (x[0]-y[0], x[1]-y[1])
ONE = (Fr(1),Fr(0))

# quotient machinery over 3 letters
letters = ["g2","g3","u"]
idx = {l:i for i,l in enumerate(letters)}
basis3 = list(itertools.product(letters, repeat=3))
bidx = {t:i for i,t in enumerate(basis3)}
Dcols = []
for a,b,c in itertools.product(letters, repeat=3):
    v = sp.zeros(27,1)
    for t in [(a,b,c),(b,a,c),(b,c,a)]:
        v[bidx[t],0] += 1
    Dcols.append(v)
DM = sp.Matrix.hstack(*Dcols)
rD = DM.rank()
print(f"letters g2,g3,u: dim V^3 = 27, rank D = {rD}, dim Q = {27-rD}")

def coord(sym):
    v = sp.zeros(27,1)
    for t,c in sym.items():
        v[bidx[t],0] += sp.Rational(c.numerator, c.denominator)
    return v
def S_Li3v(vz, vomz):
    out = {}
    for a,ca in vomz.items():
        for b,cb in vz.items():
            for c,cc in vz.items():
                out[(a,b,c)] = out.get((a,b,c), Fr(0)) - ca*cb*cc
    return out

cols, zs = [], []
for z in sorted(CR12):
    vz = F3.vec(z); vomz = F3.vec(qsub(ONE,z))
    cols.append(coord(S_Li3v(vz,vomz))); zs.append(z)
M = sp.Matrix.hstack(*cols)
r_all = sp.Matrix.hstack(M, DM).rank()
print(f"#CR = {len(zs)}: rank span S(Li3(CR)) mod D = {r_all - rD}")

# the 8 DOD5 points
z8 = [(2,-1),(Fr(-1,2),Fr(1,2)),(0,Fr(1,3)),(7,-4),(1,Fr(-1,2)),(-1,Fr(2,3)),(Fr(1,2),Fr(-1,6)),(Fr(1,3),0)]
z8 = [tuple(Fr(t) for t in z) for z in z8]
names = ["2-s3","(s3-1)/2","s3/3","7-4s3","(2-s3)/2","(2s3-3)/3","(3-s3)/6","1/3"]
allin = True
for nm,z in zip(names,z8):
    vz = F3.vec(z); vomz = F3.vec(qsub(ONE,z))
    t = coord(S_Li3v(vz,vomz))
    inspan = sp.Matrix.hstack(M, DM, t).rank() == r_all
    allin &= inspan
    print(f"  S(Li3({nm})) in span(S(Li3(CR)))+D: {inspan}")
print("ALL 8 DOD5 Li3 symbols in CR-span mod products:", allin)
