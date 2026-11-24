"""cwred.py -- exact IBP reducer for the corner family in w, generic field F.
PF = [realdict {pole: coeff}, ppos [e1..em] at +theta, pneg at -theta, const],
theta^2 = -r.  Real simple poles in {0, +-1, +-r}.  Fold F(R) = r R(r/w)/w^2.
Mirrors WTS fastred.py (v -> w, i -> theta, v^2+1 -> w^2+r, {0,+-a,+-b} -> {0,+-1,+-r}, v->1/v fold -> w->r/w).
"""
from fractions import Fraction as Fr
from math import comb

class Reducer:
    def __init__(self, F, r):
        self.F = F
        self.r = F.c(r)
        # theta: element with theta^2 = -r.  For Quad fields: theta = (0,1).  For QSI: theta = i*b.
        if hasattr(F, 'modulus'):      # Quad: theta = 0 + 1*theta
            self.th = F(0, 1)
            assert (self.th*self.th + self.r).is_zero()
        else:                           # QSI squircle: theta = i b, b = sqrt2-1; r = 3-2sqrt2
            import sympy
            self.th = F(0, -1, 1, 1)    # i(sqrt2-1) = -i + i sqrt2 -> p2 = -1+1? see below
            # i*b = i*(sqrt2-1) = -i + sqrt2 i -> p2 = -1, p3 = 1
            self.th = F(0, 0, -1, 1)
            assert (self.th*self.th + self.r).is_zero(), (self.th*self.th)
        self.ZERO = F(0); self.ONE = F(1)
        r_ = self.r
        self.roots = [F(0), F(1), F(-1), r_, -r_]
        self.DL1 = {r_: F(1), F(-1): F(1), F(0): F(-1)}    # d lt/dw
        self.DL2 = {F(1): F(1), -r_: F(1), F(0): F(-1)}    # d lmt/dw

    # ---- PF basics ----
    def pf_zero(self): return [{}, [], [], self.ZERO]
    def pf_copy(self, P): return [dict(P[0]), list(P[1]), list(P[2]), P[3]]
    def pf_add(self, P, Q):
        R = self.pf_copy(P)
        Z = self.ZERO
        for r, c in Q[0].items():
            R[0][r] = R[0].get(r, Z) + c
            if R[0][r].is_zero(): del R[0][r]
        def addlist(L1, L2):
            n = max(len(L1), len(L2)); out = []
            for t in range(n):
                x = L1[t] if t < len(L1) else Z
                y = L2[t] if t < len(L2) else Z
                out.append(x + y)
            while out and out[-1].is_zero(): out.pop()
            return out
        R[1] = addlist(R[1], Q[1]); R[2] = addlist(R[2], Q[2])
        R[3] = R[3] + Q[3]
        return R
    def pf_neg(self, P):
        return [{r: -c for r, c in P[0].items()}, [-e for e in P[1]], [-f for f in P[2]], -P[3]]
    def pf_sub(self, P, Q): return self.pf_add(P, self.pf_neg(Q))
    def pf_scale(self, P, q):
        q = self.F.c(q)
        if q.is_zero(): return self.pf_zero()
        return [{r: c*q for r, c in P[0].items()}, [e*q for e in P[1]], [f*q for f in P[2]], P[3]*q]
    def pf_check_reality(self, P, ctx=""):
        for r, c in P[0].items():
            assert r.is_real() and c.is_real(), (ctx, "real", r, c)
        Z = self.ZERO
        for e, f in zip(P[1] + [Z]*max(0, len(P[2])-len(P[1])),
                        P[2] + [Z]*max(0, len(P[1])-len(P[2]))):
            assert e == f.conj(), (ctx, "conj", e, f)
        assert P[3].is_real(), (ctx, "const", P[3])

    # ---- fold: R -> r R(r/w)/w^2 ----
    def pf_fold(self, P):
        F = self.F; rpar = self.r; th = self.th; Z = self.ZERO
        real, ppos, pneg, const = P
        assert const.is_zero(), "fold of nonzero constant"
        R = self.pf_zero()
        for rho, c in real.items():
            if rho.is_zero():
                R[0][rho] = R[0].get(rho, Z) + c
            else:
                R[0][F(0)] = R[0].get(F(0), Z) + c
                rinv = rpar / rho
                R[0][rinv] = R[0].get(rinv, Z) - c
        def fold_side(coeffs, sign):
            # pole at sign*theta (given as w - sign th), maps to w + sign th; unit u = -sign th:
            # r - sign th w = u (w + sign th)
            sth = th if sign > 0 else -th
            u = -sth
            out_real0 = Z
            out = []
            def add_at(order, val):
                while len(out) < order: out.append(Z)
                out[order-1] = out[order-1] + val
            for j, e in enumerate(coeffs, start=1):
                if e.is_zero(): continue
                if j == 1:
                    c0 = e * rpar / u / sth     # = e
                    out_real0 = out_real0 + c0
                    add_at(1, -c0)
                else:
                    for s in range(0, j-1):
                        coef = e * rpar * Fr(comb(j-2, s)) * ((-sth) ** (j-2-s)) / (u ** j)
                        add_at(j - s, coef)
            return out_real0, out
        r0a, newneg = fold_side(ppos, +1)
        r0b, newpos = fold_side(pneg, -1)
        if not (r0a + r0b).is_zero():
            R[0][F(0)] = R[0].get(F(0), Z) + r0a + r0b
        R[1] = newpos; R[2] = newneg
        for L in (R[1], R[2]):
            while L and L[-1].is_zero(): L.pop()
        for rr, c in list(R[0].items()):
            if c.is_zero(): del R[0][rr]
        return R

    def pf_split(self, P):
        real, ppos, pneg, const = P
        simple = self.pf_zero(); higher = self.pf_zero()
        simple[0] = dict(real)
        if ppos: simple[1] = [ppos[0]]
        if pneg: simple[2] = [pneg[0]]
        Z = self.ZERO
        higher[1] = [Z] + ppos[1:] if len(ppos) >= 2 else []
        higher[2] = [Z] + pneg[1:] if len(pneg) >= 2 else []
        for H in (higher[1], higher[2]):
            while H and H[-1].is_zero(): H.pop()
        return simple, higher, const

    def pf_integrate_higher(self, H):
        def integ(coeffs):
            out = []
            for j, e in enumerate(coeffs, start=1):
                if j == 1:
                    assert e.is_zero(), "log term in higher antiderivative"
                    continue
                out.append(e.scale(Fr(-1, j-1)))
            while out and out[-1].is_zero(): out.pop()
            return out
        G = self.pf_zero()
        G[1] = integ(H[1]); G[2] = integ(H[2])
        return G

    def pf_eval_at(self, G, x):
        x = self.F.c(x)
        tot = self.ZERO; th = self.th
        for j, e in enumerate(G[1], start=1):
            tot = tot + e / ((x - th) ** j)
        for j, f in enumerate(G[2], start=1):
            tot = tot + f / ((x + th) ** j)
        tot = tot + G[3]
        for r, c in G[0].items():
            tot = tot + c / (x - r)
        return tot

    def pf_mul_DL(self, Gm, DL):
        realG, ppos, pneg, const = Gm
        Z = self.ZERO; th = self.th
        R = self.pf_zero()
        for rho, d in DL.items():
            if not const.is_zero():
                R[0][rho] = R[0].get(rho, Z) + const * d
            for coeffs, sgn in ((ppos, +1), (pneg, -1)):
                si = th if sgn > 0 else -th
                den = rho - si
                for j, e in enumerate(coeffs, start=1):
                    if e.is_zero(): continue
                    ed = e * d
                    R[0][rho] = R[0].get(rho, Z) + ed / (den ** j)
                    for u in range(1, j+1):
                        val = -ed / (den ** (j - u + 1))
                        L = R[1] if sgn == +1 else R[2]
                        while len(L) < u: L.append(Z)
                        L[u-1] = L[u-1] + val
        for rho, c in list(R[0].items()):
            if c.is_zero(): del R[0][rho]
        for L in (R[1], R[2]):
            while L and L[-1].is_zero(): L.pop()
        return R

    # ---- from numerator/denominator polynomials over F ----
    def pf_from_poly(self, N, D):
        F = self.F; Z = self.ZERO; ONE = self.ONE; th = self.th; rpar = self.r
        def peval(poly, x):
            tot = Z
            for c in reversed(poly): tot = tot * x + c
            return tot
        def pdeflate(poly, x):
            n = len(poly) - 1
            q = [Z]*n
            carry = poly[n]; q[n-1] = carry
            for t in range(n-1, 0, -1):
                carry = poly[t] + carry * x; q[t-1] = carry
            rem = poly[0] + carry * x
            assert rem.is_zero(), ("deflation remainder", x, rem)
            return q
        def pmul(P1, P2):
            out = [Z]*(len(P1)+len(P2)-1)
            for a, ca in enumerate(P1):
                for b, cb in enumerate(P2):
                    out[a+b] = out[a+b] + ca*cb
            return out
        def ptrimz(P):
            while P and P[-1].is_zero(): P.pop()
            return P
        N = ptrimz(list(N)); D = ptrimz(list(D))
        assert len(D) >= 1
        P = self.pf_zero()
        if len(N) >= len(D):
            q = [Z]*(len(N)-len(D)+1)
            rem = list(N)
            for t in range(len(N)-len(D), -1, -1):
                c = rem[len(D)-1+t] / D[-1]
                q[t] = c
                for s in range(len(D)):
                    rem[s+t] = rem[s+t] - c*D[s]
            while rem and rem[-1].is_zero(): rem.pop()
            assert len(q) <= 1, ("nonconstant polynomial part", q)
            if q: P[3] = q[0]
            N = rem
        Drem = list(D)
        roots_found = []
        for rho in self.roots:
            while len(Drem) >= 2 and peval(Drem, rho).is_zero():
                Drem2 = pdeflate(Drem, rho)
                if peval(Drem2, rho).is_zero():
                    raise ValueError(f"higher real pole at {rho}")
                res = peval(N, rho) / peval(pdeflate(list(D), rho), rho)
                P[0][rho] = P[0].get(rho, Z) + res
                Drem = Drem2
                roots_found.append(rho)
        for rr, c in list(P[0].items()):
            if c.is_zero(): del P[0][rr]
        m = (len(Drem) - 1)//2
        assert 2*m + 1 == len(Drem), ("odd leftover denominator", Drem)
        lead = Drem[-1]
        pw = [ONE]
        for _ in range(m):
            new = [Z]*(len(pw)+2)
            for t, c in enumerate(pw):
                new[t] = new[t] + c * rpar
                new[t+2] = new[t+2] + c
            pw = new
        for a1, b1 in zip(pw, Drem):
            assert (a1*lead - b1).is_zero(), ("denominator not lead*(w^2+r)^m", Drem)
        H = [lead]
        for _ in range(m):
            H = pmul(H, [2*th, ONE])            # (z + 2 th)^m * lead
        for rho in roots_found:
            H = pmul(H, [th - rho, ONE])        # (z + (th - rho))
        Nw = []
        for s in range(m):
            tot = Z
            for t, c in enumerate(N):
                if t >= s:
                    tot = tot + c * Fr(comb(t, s)) * (th ** (t-s))
            Nw.append(tot)
        q = []
        for s in range(m):
            acc = Nw[s]
            for t in range(1, s+1):
                acc = acc - (H[t] if t < len(H) else Z) * q[s-t]
            q.append(acc / H[0])
        ppos = [q[m-j] for j in range(1, m+1)]
        while ppos and ppos[-1].is_zero(): ppos.pop()
        P[1] = ppos
        P[2] = [e.conj() for e in ppos]
        return P

    # ---- the reduction ----
    def fold_key(self, key, R):
        i, j = key
        if i >= j:
            return key, R
        return (j, i), self.pf_fold(R)

    def reduce_k(self, integrand_pf):
        work = {}
        for key, R in integrand_pf.items():
            kk, RR = self.fold_key(key, R)
            work[kk] = self.pf_add(work.get(kk, self.pf_zero()), RR)
        final = {}
        kmax = max(sum(kk) for kk in work)
        for deg in range(kmax, -1, -1):
            for key in sorted([kk for kk in work if sum(kk) == deg], reverse=True):
                R = work.pop(key)
                if deg == 0:
                    final[key] = R
                    continue
                simple, higher, const = self.pf_split(R)
                assert const.is_zero(), f"polynomial part at {key}"
                if simple[0] or simple[1] or simple[2]:
                    final[key] = self.pf_add(final.get(key, self.pf_zero()), simple)
                if not higher[1] and not higher[2]:
                    continue
                G = self.pf_integrate_higher(higher)
                C = self.pf_eval_at(G, self.r)
                Gm = [G[0], G[1], G[2], -C]
                i, j = key
                if i >= 1:
                    kk, RR = self.fold_key((i-1, j), self.pf_scale(self.pf_mul_DL(Gm, self.DL1), Fr(-i)))
                    work[kk] = self.pf_add(work.get(kk, self.pf_zero()), RR)
                if j >= 1:
                    kk, RR = self.fold_key((i, j-1), self.pf_scale(self.pf_mul_DL(Gm, self.DL2), Fr(-j)))
                    work[kk] = self.pf_add(work.get(kk, self.pf_zero()), RR)
        for key in final:
            self.pf_check_reality(final[key], ctx=str(key))
        return final
