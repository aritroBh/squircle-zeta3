# step_24: fast exact quadratic-field vecfun (Fraction-pair arithmetic)
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
        """letters: (name, (a,b), norm_int). Every z must be +- prod letters^{a_i}, a_i in Z."""
        self.D = Fr(D); self.letters = letters
    def elem(self, v):
        r = (Fr(1), Fr(0))
        for nm, e in v.items():
            lt = dict((n,(ab,nm_)) for n,ab,nm_ in self.letters)[nm]
            r = qmul(r, qpow(lt, int(e), self.D), self.D)
        return r
    def vec(self, z):
        """z: (a,b) exact a+b sqrt(D). Returns dict letter->int exp (exact, verified)."""
        a,b = Fr(z[0]), Fr(z[1])
        N = a*a - self.D*b*b
        absN = abs(N)
        exps = {}
        rem_norm = absN
        # norm-driven exponents (order: letters with |norm|>1 first)
        for nm, ab, nrm in self.letters:
            if abs(nrm) == 1: continue
            base = abs(nrm)
            # count exponent of prime(s) of base in rem_norm
            e = 0
            while rem_norm % base == 0 and rem_norm > 1:
                rem_norm //= base; e += 1
            if nrm < 0 and abs(nrm) == 2 and self.D == 2:
                pass
            exps[nm] = e
        if rem_norm != 1:
            raise ValueError(f"norm leftover {rem_norm} for {z}")
        # unit exponent via logs (fast, then exact verify)
        mp.mp.dps = 50
        target = abs(mp.mpf(int(a)) + mp.sqrt(int(self.D))*mp.mpf(int(b))) if b else abs(mp.mpf(int(a)))
        lt = mp.log(target)
        prod_log = mp.mpf(0)
        units = []
        for nm, ab, nrm in self.letters:
            e = exps.get(nm, 0)
            if e:
                val = abs(mp.mpf(int(ab[0])) + mp.sqrt(int(self.D))*mp.mpf(int(ab[1]))) if ab[1] else abs(mp.mpf(int(ab[0])))
                prod_log += e*mp.log(val)
        for nm, ab, nrm in self.letters:
            if abs(nrm) == 1: units.append(nm)
        resid = lt - prod_log
        if units:
            un = units[0]
            uab = dict((n,(ab2,nm_)) for n,ab2,nm_ in self.letters)[un]
            uval = abs(mp.mpf(int(uab[0])) + mp.sqrt(int(self.D))*mp.mpf(int(uab[1])))
            c = resid/mp.log(uval)
            ci = int(mp.nint(c))
            assert abs(c-ci) < mp.mpf('1e-30'), (z, c)
            exps[un] = exps.get(un, 0) + ci
        # EXACT verify
        got = self.elem(exps)
        if got != (a,b) and got != (-a,-b):
            raise ValueError(f"verify fail {z} -> {exps} gives {got}")
        return {k:Fr(v) for k,v in exps.items() if v}

F2 = QF(2, [("r2",(Fr(0),Fr(1)),-2), ("u",(Fr(1),Fr(1)),-1)])
F3 = QF(3, [("g2",(Fr(1),Fr(1)),-2), ("g3",(Fr(0),Fr(1)),-3), ("u",(Fr(2),Fr(1)),1)])
F5 = QF(5, [("2",(Fr(2),Fr(0)),4), ("s5",(Fr(0),Fr(1)),-5), ("phi",(Fr(1,2),Fr(1,2)),-1)])

# test on the 12-gon points and p_j
s3tests = [( (2,-1)), ( (Fr(-1,2),Fr(1,2)) ), ( (0,Fr(1,3)) ), ( (7,-4) ), ( (1,Fr(-1,2)) ),
           ( (-1,Fr(2,3)) ), ( (Fr(1,2),Fr(-1,6)) ), ( (Fr(1,3),0) ),
           ((3,1)), ((4,2)), ((-3,-2)), ((-1,Fr(-1,2))) ]
import sympy as sp
def pairname(z):
    return f"{z[0]}+{z[1]}*sqrt3"
for z in s3tests:
    try:
        print(pairname(z), "->", F3.vec(z))
    except Exception as ex:
        print(pairname(z), "FAIL", ex)
