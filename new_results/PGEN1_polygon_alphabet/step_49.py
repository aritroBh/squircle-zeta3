# step_25: fast exact quadratic-field vecfun, fixed
from fractions import Fraction as Fr
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

class QF:
    def __init__(self, D, letters):
        """letters: list of (name, (a,b), norm)."""
        self.D = Fr(D); self.letters = letters
    def elem(self, v):
        r = (Fr(1), Fr(0))
        L = {n: ab for n, ab, _ in self.letters}
        for nm, e in v.items():
            r = qmul(r, qpow(L[nm], int(e), self.D), self.D)
        return r
    def vec(self, z):
        a,b = Fr(z[0]), Fr(z[1])
        N = a*a - self.D*b*b
        rem = abs(N)
        exps = {}
        for nm, ab, nrm in self.letters:
            if abs(nrm) == 1: continue
            base = abs(nrm); e = 0
            if rem >= 1:
                while rem % base == 0: rem //= base; e += 1
            else:
                while rem < 1: rem *= base; e -= 1
                # after this rem might be >1 and divisible further: handle by re-loop both ways
            exps[nm] = e
        # second pass for mixed (e.g. rem grew past 1)
        changed = True
        while changed:
            changed = False
            for nm, ab, nrm in self.letters:
                if abs(nrm) == 1: continue
                base = abs(nrm)
                while rem % base == 0 and rem != 1:
                    rem //= base; exps[nm] += 1; changed = True
                while rem < 1:
                    rem *= base; exps[nm] -= 1; changed = True
        if rem != 1: raise ValueError(f"norm leftover {rem} for {z}")
        # unit exponent via logs, then exact verify
        mp.mp.dps = 40
        sq = mp.sqrt(int(self.D))
        def aval(ab):
            return abs(mp.mpf(int(ab[0])) + sq*mp.mpf(int(ab[1]))) if ab[1] else abs(mp.mpf(int(ab[0])))
        lt = mp.log(aval((a,b)))
        prod_log = mp.mpf(0)
        for nm, ab, nrm in self.letters:
            e = exps.get(nm, 0)
            if e and abs(nrm) != 1: prod_log += e*mp.log(aval(ab))
        units = [ (nm,ab) for nm,ab,nrm in self.letters if abs(nrm)==1 ]
        resid = lt - prod_log
        if units:
            nm, ab = units[0]
            c = resid/mp.log(aval(ab))
            ci = int(mp.nint(c))
            if abs(c-ci) > mp.mpf('1e-25'): raise ValueError(f"unit exp not integer: {z} {c}")
            exps[nm] = exps.get(nm,0)+ci
        if abs(resid) > 1e-30 and not units:
            raise ValueError("residual log nonzero")
        got = self.elem(exps)
        if got != (a,b) and got != (-a,-b):
            raise ValueError(f"verify fail {z} -> {exps} gives {got}")
        return {k:Fr(v) for k,v in exps.items() if v}

F3 = QF(3, [("g2",(Fr(1),Fr(1)),-2), ("g3",(Fr(0),Fr(1)),-3), ("u",(Fr(2),Fr(1)),1)])
F2 = QF(2, [("r2",(Fr(0),Fr(1)),-2), ("u",(Fr(1),Fr(1)),-1)])
F5 = QF(5, [("2",(Fr(2),Fr(0)),4), ("s5",(Fr(0),Fr(1)),-5), ("phi",(Fr(1,2),Fr(1,2)),-1)])

tests3 = [ (2,-1), (Fr(-1,2),Fr(1,2)), (0,Fr(1,3)), (7,-4), (1,Fr(-1,2)),
           (-1,Fr(2,3)), (Fr(1,2),Fr(-1,6)), (Fr(1,3),0),
           (3,1), (4,2), (-3,-2), (-1,Fr(-1,2)), (Fr(-1,2),Fr(-1,2)), (Fr(3,2),Fr(-1,2)) ]
for z in tests3:
    print(z, "->", F3.vec(z))
print("F3 OK")
