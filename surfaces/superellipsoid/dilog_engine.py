# dilog_engine.py -- exact reduction of real dilogarithms at S-units of Q(sqrt2), S={2,3},
# using only classical functional equations (reflection, Landen, inversion, duplication,
# Abel's five-term relation on {x<1,y<1,x+y<1}).  Every relation used is an explicit instance;
# the linear algebra is exact (Fractions).  A certificate (list of instances with multipliers)
# can be extracted for any reduction.
from fractions import Fraction as Fr
from qs2 import (Q2, W2, Li2, PI2, R2, A, B, factor, rel_reflection, rel_landen, rel_inversion,
                 rel_dup, rel_abel, canon_zero_li2)
import itertools, math

def pool(maxi=6, maxj=2, maxk=6):
    pts = set()
    for i in range(-maxi, maxi+1):
        for j in range(-maxj, maxj+1):
            for k in range(-maxk, maxk+1):
                for sg in (1, -1):
                    x = Q2(sg)*(R2**i)*(Q2(3)**j)*(A**k)
                    fx = float(x)
                    if not (fx < 1) or x == Q2(0): continue
                    if factor(Q2(1) - x) is None: continue
                    pts.add(x)
    return pts

def relations(P, abel=True):
    rels = []   # (name, W2)
    for x in P:
        fx = float(x)
        if 0 < fx < 1 and (Q2(1)-x) in P:
            rels.append((("R", x), rel_reflection(x)))
        y = x/(x - Q2(1))
        if y in P and fx != 0:
            rels.append((("La", x), rel_landen(x)))
        if fx < 0 and (Q2(1)/x) in P:
            rels.append((("I", x), rel_inversion(x)))
        if abs(fx) <= 1 and (-x) in P and (x*x) in P:
            rels.append((("D", x), rel_dup(x)))
    if abel:
        L = sorted(P, key=float)
        for x, y in itertools.combinations_with_replacement(L, 2):
            if not (float(x) + float(y) < 1): continue
            if factor(Q2(1) - x - y) is None: continue
            args = [x/(Q2(1)-y), y/(Q2(1)-x), x*y/((Q2(1)-x)*(Q2(1)-y))]
            if all((t in P) or t == Q2(0) for t in args):
                rels.append((("A", x, y), rel_abel(x, y)))
    return [(n, canon_zero_li2(w)) for n, w in rels]

def colkey(k):
    return k

class Reducer:
    def __init__(self, P, basis, abel=True, verbose=True):
        self.P = P
        self.basis = [Q2(1) * b for b in basis]
        rels = relations(P, abel)
        if verbose: print("pool size", len(P), " relations", len(rels))
        # column order: non-basis Li2 first (eliminate), basis Li2 later, constants never pivot
        bset = set(self.basis)
        def prio(key):
            if key[0] == 'Li2':
                return (0 if key[1] not in bset else 1, float(key[1]))
            return (2, str(key))
        self.prio = prio
        self.rows = {}     # pivot key -> (row dict of Fractions, certificate dict)
        for name, w in rels:
            row = {}
            for k, v in w.d.items():
                if v.q != 0: raise ValueError("relation with irrational coefficient?")
                row[k] = v.p
            self._insert(row, {name: Fr(1)})
        if verbose: print("echelon rows", len(self.rows))
    def _reduce_row(self, row, cert):
        row = dict(row); cert = dict(cert)
        while True:
            keys = [k for k in row if k[0] == 'Li2' and k in self.rows]
            if not keys: return row, cert
            k = min(keys, key=self.prio)
            prow, pcert = self.rows[k]
            f = row[k]/prow[k]
            for kk, vv in prow.items():
                nv = row.get(kk, Fr(0)) - f*vv
                if nv == 0: row.pop(kk, None)
                else: row[kk] = nv
            for kk, vv in pcert.items():
                nv = cert.get(kk, Fr(0)) - f*vv
                if nv == 0: cert.pop(kk, None)
                else: cert[kk] = nv
    def _insert(self, row, cert):
        row, cert = self._reduce_row(row, cert)
        li = [k for k in row if k[0] == 'Li2']
        if not li: return
        piv = min(li, key=self.prio)
        self.rows[piv] = (row, cert)
    def reduce(self, w):
        """reduce a W2 (Q2 coefficients) -> (W2 reduced, certificate for rational part, for sqrt2 part)"""
        out = W2(); certs = []
        for part in (0, 1):
            row = {}
            for k, v in w.d.items():
                c = v.p if part == 0 else v.q
                if c != 0: row[k] = Fr(c)
            red, cert = self._reduce_row(row, {})
            for k, c in red.items():
                out.add(k, Q2(c) if part == 0 else Q2(0, c))
            certs.append(cert)
        return out, certs

if __name__ == "__main__":
    import mpmath as mp
    mp.mp.dps = 50
    P = pool()
    red = Reducer(P, basis=[B, -(B**4), Q2(1, 0)/3])
    # tests: note's Lemma 3 identities, e.g. Li2(-b) and Li2(b^2), Li2(2-sqrt2), Li2(-b^2), Li2((2+sqrt2)/4)
    for x in [-B, B*B, Q2(2, -1), -B*B, Q2(Fr(1, 2), Fr(1, 4)), Q2(-1, 0)/2, Q2(1, 0)/4, Q2(2, 0)/3, R2*B/3, A/3, -(B**4)]:
        if x not in P: print("not in pool", x); continue
        r, c = red.reduce(Li2(x))
        print("Li2(%r) = %r   [check %s]" % (x, r, mp.nstr(mp.polylog(2, x.mp()) - r.num(), 3)))
