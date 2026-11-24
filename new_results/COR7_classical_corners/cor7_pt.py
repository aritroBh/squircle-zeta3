"""cor7_pt.py -- exact letter-vector machinery for quadratic corners (class number 1),
with guided-then-exactly-certified exponent vectors (no log-lattice inversion).

F = Q(sqrt D), D squarefree, h(D) = 1. For the letter set we build generators:
  eps   : fundamental unit (N = +-1, smallest |log| > 0)
  for each rational prime l dividing any letter norm:
      inert    -> generator l itself (N = l^2)
      ramified -> element g_l with N(g_l) = +-l
      split    -> elements l (N = l^2) and g_l (N(g_l) = +-l), one prime above l
vec(x): x assumed S-unit: norm valuations constrain prime exponents to a small range;
unit exponent from one embedding (rounded); every result EXACTLY verified in NF arithmetic,
else None (=> the letter is not in the generated group).
"""
from fractions import Fraction as Fr
import itertools, sys
import mpmath as mp
sys.path.insert(0, '<project>/research/08_theory/COR6_true_point_symbols')
from cor6_nf import NF
from cor6_sym import sym_J, sscale
from cor6_ft import ft_add, ft_scale

def _sqf(D):
    f = {}
    d, n = 2, D
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1 if d == 2 else 2
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def factorint(n):
    return _sqf(n)

class QuadField:
    def __init__(self, D, verbose=False):
        self.D = D
        self.F = NF([-D, 0, 1])
        mp.mp.dps = 80
        self.rD = mp.sqrt(D)
        self.eps = self._find_unit()
        self.gen_primes = {}   # l -> list of (element, |N| power of l: 1 or 2)
        self.mode = {}
        if verbose:
            print(f"Q(sqrt{D}): eps = {self.eps}")
    def _find_unit(self, bound=200):
        best = None
        for e in (1, 2):
            if e == 2 and self.D % 4 != 1:
                continue
            for a in range(-bound, bound+1):
                for b in range(1, bound+1):
                    if abs(a*a - self.D*b*b) == e*e:
                        v = abs(a + b*float(self.rD))/e
                        if v > 0:
                            lg = abs(mp.log(v))
                            if lg > 1e-9 and (best is None or lg < best[0]):
                                best = (lg, (Fr(a, e), Fr(b, e)))
        return best[1] if best else None
    def _find_norm_l(self, l, bound=200):
        """element with |N| = l (smallest log), or None."""
        best = None
        for e in (1, 2):
            if e == 2 and self.D % 4 != 1:
                continue
            for a in range(-bound, bound+1):
                for b in range(1, bound+1):
                    if abs(a*a - self.D*b*b) == l*e*e:
                        v = abs(a + b*float(self.rD))/e
                        if v > 0:
                            lg = abs(mp.log(v))
                            if best is None or lg < best[0]:
                                best = (lg, (Fr(a, e), Fr(b, e)))
        return best[1] if best else None
    def prime_status(self, l):
        if l in self.mode:
            return self.mode[l]
        r = sp_legendre(self.D, l)
        st = 'ram' if r == 0 else ('split' if r == 1 else 'inert')
        self.mode[l] = st
        return st
    def generators_for(self, primes):
        """build generator list [(elem, kind, l)] covering the given rational primes."""
        gens = [('eps', self.eps, 1)]
        out = []
        for l in sorted(primes):
            st = self.prime_status(l)
            if st == 'inert':
                out.append((self.F.from_int(l), 'inert', l))
            else:
                g = self._find_norm_l(l)
                if g is None:
                    return None
                if st == 'ram':
                    out.append((g, 'ram', l))
                else:
                    out.append((self.F.from_int(l), 'splitQ', l))
                    out.append((g, 'splitp', l))
        return gens + out

def sp_legendre(D, l):
    """splitting of rational prime l in Q(sqrt D), D squarefree: 'split'/'inert'/'ram' encoded 1/-1/0."""
    if l == 2:
        if D % 2 == 0: return 0
        r = D % 8
        if r == 1: return 1
        if r == 5: return -1
        return 0   # D == 3 mod 4: disc 4D, 2 ramifies
    if D % l == 0: return 0
    return 1 if pow(D % l, (l-1)//2, l) == 1 else -1

class QuadPoint:
    def __init__(self, QF, label, c0, primes):
        self.QF = QF
        self.F = QF.F
        self.label = label
        self.c0 = c0
        self.gens = QF.generators_for(primes)
        if self.gens is None:
            raise ValueError("no principal generator for some prime")
        self.ginfo = []
        for kind, g, l in self.gens:
            pass
        # normalize to (elem, kind, l)
        self.gensE = []
        for item in self.gens:
            if item[0] == 'eps':
                self.gensE.append((item[1], 'eps', 1))
            else:
                self.gensE.append(item)
        # norms
        self.gnorm = []
        for el, kind, l in self.gensE:
            self.gnorm.append(self.F.norm(el))
        self.cache = {}
        zero = self.F.zero(); one = self.F.one()
        self.L = {0: zero, 1: one, -1: self.F.neg(one)}
    def vec(self, el):
        if el in self.cache:
            return self.cache[el]
        F = self.F
        N = F.norm(el)
        if N == 0:
            raise ValueError("zero")
        vN = {}
        n, d = abs(N).numerator, abs(N).denominator
        for p, e in (factorint(n) if n > 1 else {}).items(): vN[p] = vN.get(p, 0) + e
        for p, e in (factorint(d) if d > 1 else {}).items(): vN[p] = vN.get(p, 0) - e
        gens = self.gensE
        mp.mp.dps = 100
        def logv(g, i=0):
            r = -self.QF.rD if i == 0 else self.QF.rD
            s = mp.mpf(0); t = mp.mpf(1)
            for cf in g:
                s += mp.mpf(cf.numerator)/cf.denominator * t
                t *= r
            return mp.log(abs(s)) if s != 0 else None
        lx = logv(el, 0)
        logs = [logv(g[0], 0) for g in gens]
        eps_i = next(i for i in range(len(gens)) if gens[i][1] == 'eps')
        # enumerate splitp exponents; others determined by norm valuations
        idx_splitp = [i for i in range(len(gens)) if gens[i][1] == 'splitp']
        ranges = []
        for i in idx_splitp:
            l_ = gens[i][2]
            v = vN.get(l_, 0)
            ranges.append(range(v - 16, v + 17))
        for combo in itertools.product(*ranges) if ranges else [()]:
            e = [0]*len(gens)
            for j, i in enumerate(idx_splitp):
                e[i] = combo[j]
            ok = True
            for i in range(len(gens)):
                el_, k_, l_ = gens[i]
                if k_ == 'eps' or k_ == 'splitp':
                    continue
                v = vN.get(l_, 0)
                for j in idx_splitp:
                    if gens[j][2] == l_:
                        v -= e[j]
                if k_ == 'ram':
                    e[i] = v
                elif k_ == 'inert':
                    if v % 2: ok = False; break
                    e[i] = v//2
                elif k_ == 'splitQ':
                    if v % 2: ok = False; break
                    e[i] = v//2
            if not ok:
                continue
            rem = lx
            for i in range(len(gens)):
                if i != eps_i and e[i]:
                    rem -= e[i]*logs[i]
            e[eps_i] = int(mp.nint(rem / logs[eps_i]))
            v = F.one()
            for i in range(len(gens)):
                if e[i]:
                    v = F.mul(v, F.pow(gens[i][0], e[i]))
            if F.eq(v, el) or F.eq(v, F.neg(el)):
                out = tuple(e)
                self.cache[el] = out
                return out
        self.cache[el] = None
        return None
    def letter_span_rank(self):
        F, one = self.F, self.F.one()
        labs = [("2", F.from_int(2)), ("c", self.c0),
                ("1-c", F.sub(one, self.c0)), ("1+c", F.add(one, self.c0))]
        vecs = []
        for nm, el in labs:
            v = self.vec(el)
            if v is None:
                return None, labs, None
            vecs.append((nm, v))
        M = [list(map(Fr, v)) for _, v in vecs]
        R = 0
        for c in range(len(M[0])):
            p_ = next((r for r in range(R, len(M)) if M[r][c]), None)
            if p_ is None: continue
            M[R], M[p_] = M[p_], M[R]
            q = M[R][c]; M[R] = [x/q for x in M[R]]
            for r in range(len(M)):
                if r != R and M[r][c]:
                    M[r] = [a - b*M[R][c] for a, b in zip(M[r], M[R])]
            R += 1
        return R, labs, vecs
    def sym_H(self, word):
        wordF = tuple(self.L[a] for a in word)
        s = -1 if sum(1 for a in word if a == 1) % 2 else 1
        return sscale(s, sym_J(wordF, self.c0, self.vec, a0=self.F.zero(),
                               sub=self.F.sub, iszero=lambda a: all(x == 0 for x in a)))
    def omega_symbol(self, W4TERMS):
        SOm = {}
        for e, w, val in W4TERMS:
            ce = self.F.pow(self.c0, e)
            SOm = ft_add(SOm, ft_scale(self.F, self.F.scale(val, ce), self.sym_H(w)))
        return SOm
