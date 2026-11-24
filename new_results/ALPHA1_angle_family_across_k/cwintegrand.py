"""cwintegrand.py -- exact corner-family integrand in the w variable, generic field.

Psi_eps = (sqrt(t^{2e-2} + (1-t)^{2e-2} - 2 ca (t(1-t))^{e-1}) - t^{e-1} - (1-t)^{e-1})/2, ca = cos alpha.
[eps^k] Psi = sum_{m=0}^k R~_m(t) x^m y^{k-m}, x = lt, y = lmt:
  sqrt part: R~_m = (1/2) q_w T_m(t) / (t(1-t) q2a^k),  q2a = 1 - 4 c^2 t(1-t) = t^2+(1-t)^2 - 2 ca t(1-t),
  T_m(t) = sum_{j>=1} binom(1/2,j) sum_{s1+s2+s3=j} multinom * inner * (1-t)^{2s1} t^{2s2} (t(1-t))^{s3} (-2ca)^{s3} q2a^{k-j},
  inner = sum_{m1+m2+m3=k} w(s1,m1) w(s2,m2) w'(s3,m3) binom(m3, m-m1),
  w(s,m) = s! S(m,s) 2^m/m!,  w'(s,m) = s! S(m,s)/m!.
  A part (m=k): -1/(2 k! t);  B part (m=0): -1/(2 k! (1-t)).

w-map: t = (w-r)(w+1)/(2(1-r) w), r = (1-c)/(1+c), c = cos(alpha/2); q2a = (w^2+r)^2/((1+r)^2 w^2);
R_m(w) (with dw) = C_k * T~_m(w) / [w (w^2+r)^{2k-2} (w-r)(w+1)(1-w)(w+r)],  C_k = 1/(2 (2c)^{2k-1}),
T~_m(w) = (2(1-r) w)^{2k} T_m(t(w)) (polynomial).
AB: R^A = -(w^2+r)/(2 k! w (w-r)(w+1)),  R^B = -(w^2+r)/(2 k! w (1-w)(w+r)).
"""
from fractions import Fraction as Fr
from math import comb, factorial

def stirling2(n, k):
    if n < 0 or k < 0: return 0
    if n == 0 and k == 0: return 1
    if k == 0 or n == 0: return 0
    prev = [0]*(k+1); prev[0] = 1
    for i in range(1, n+1):
        cur = [0]*(k+1)
        for j in range(1, min(i,k)+1):
            cur[j] = j*prev[j] + prev[j-1]
        prev = cur
    return prev[k] if k <= n else 0

_S = {}
def S2(n, k):
    if (n,k) not in _S: _S[(n,k)] = stirling2(n,k)
    return _S[(n,k)]

def w2(s, m):
    if s == 0: return Fr(1) if m == 0 else Fr(0)
    if m < s: return Fr(0)
    return Fr(factorial(s)*S2(m,s)*2**m, factorial(m))

def w1(s, m):
    if s == 0: return Fr(1) if m == 0 else Fr(0)
    if m < s: return Fr(0)
    return Fr(factorial(s)*S2(m,s), factorial(m))

_HALF = Fr(1,2)
def binom_half(j):
    r = Fr(1)
    for t in range(j): r *= (_HALF - t)
    return r / factorial(j)

# polynomials over Fractions, low-to-high
def padd(P, Q):
    n = max(len(P), len(Q)); out = [Fr(0)]*n
    for i,c in enumerate(P): out[i] += c
    for i,c in enumerate(Q): out[i] += c
    return out
def pscale(P, q): return [c*q for c in P]
def pmul(P, Q):
    if not P or not Q: return []
    out = [Fr(0)]*(len(P)+len(Q)-1)
    for a,ca in enumerate(P):
        if ca == 0: continue
        for b,cb in enumerate(Q):
            if cb: out[a+b] += ca*cb
    return out
def pshift(P, n): return [Fr(0)]*n + P
def ptrim(P):
    while P and P[-1] == 0: P.pop()
    return P
def p_pow(P, n):
    out = [Fr(1)]
    for _ in range(n): out = pmul(out, P)
    return out
def p_comp_1mt(P):
    out = [Fr(0)]
    for s, c in enumerate(P):
        if c == 0: continue
        term = p_pow([Fr(1), Fr(-1)], s)
        out = padd(out, pscale(term, c))
    return ptrim(out)

class CornerIntegrand:
    """Builds [eps^k] Psi coefficients R_m(w) as PF data (numerator/denominator polynomial pairs over F)."""
    def __init__(self, F, c, r):
        """F: field class; c, r: field elements (r = (1-c)/(1+c)). ca = 2c^2-1 must be rational."""
        self.F = F
        self.c = F.c(c); self.r = F.c(r)
        ca = self.c*self.c*2 - 1
        assert ca.is_real(), "ca must be rational for the T-builder"
        # extract rational ca
        self.ca = ca.a if hasattr(ca, 'a') else ca.p0
        self.q2a = ptrim([Fr(1), -2*(1+self.ca), 2*(1+self.ca)])   # q2alpha(t)

    def inner_sum(self, k, m, s1, s2, s3):
        tot = Fr(0)
        for m1 in range(s1, k+1):
            wA = w2(s1, m1)
            if wA == 0: continue
            for m3 in range(s3, k+1-m1):
                m2 = k - m1 - m3
                if m2 < s2: continue
                nu = m - m1
                if nu < 0 or nu > m3: continue
                wB = w2(s2, m2); wC = w1(s3, m3)
                if wB == 0 or wC == 0: continue
                tot += wA*wB*wC*comb(m3, nu)
        return tot

    def T_poly(self, k, m):
        """T_m(t) as Fraction poly."""
        out = []
        q2a = self.q2a; ca = self.ca
        for j in range(1, k+1):
            bh = binom_half(j)
            q2pow = p_pow(q2a, k-j)
            for s1 in range(j+1):
                for s3 in range(j+1-s1):
                    s2 = j - s1 - s3
                    inn = self.inner_sum(k, m, s1, s2, s3)
                    if inn == 0: continue
                    coef = bh * comb(j, s1) * comb(j-s1, s3) * inn * (-2*ca)**s3
                    if coef == 0: continue
                    # (1-t)^{2 s1} t^{2 s2} (t(1-t))^{s3} q2a^{k-j}
                    term = p_pow([Fr(1), Fr(-1)], 2*s1)
                    term = pshift(term, 2*s2)
                    term = pmul(term, pshift(p_pow([Fr(1), Fr(-1)], s3), s3))
                    term = pmul(term, q2pow)
                    out = padd(out, pscale(term, coef))
        return ptrim(out)

    # ---- w-map (polynomials over F) ----
    def fpmul(self, P, Q):
        Z = self.F(0)
        out = [Z]*(len(P)+len(Q)-1)
        for a,ca in enumerate(P):
            for b,cb in enumerate(Q):
                out[a+b] = out[a+b] + ca*cb
        while out and out[-1].is_zero(): out.pop()
        return out
    def fpscale(self, P, q): return [c*q for c in P]
    def fpadd(self, P, Q):
        Z = self.F(0)
        n = max(len(P), len(Q)); out = [Z]*n
        for i,c in enumerate(P): out[i] = out[i] + c
        for i,c in enumerate(Q): out[i] = out[i] + c
        return out
    def fp_pow(self, P, n):
        out = [self.F(1)]
        for _ in range(n): out = self.fpmul(out, P)
        return out

    def sqrt_piece_nd(self, k, T):
        """returns (num, den) over F for R_m(w)."""
        F = self.F; r = self.r
        Nt = self.fpmul([-r, F(1)], [F(1), F(1)])          # (w-r)(w+1) = w^2 + (1-r) w - r
        W = [F(0), F(1)]                                    # w
        D = len(T)-1
        # T~ = sum_d tau_d Nt^d (2(1-r) w)^{2k-d}
        Tt = [F(0)]
        fac = (F(1)-r)*2
        for d, tau in enumerate(T):
            if tau == 0: continue
            term = self.fp_pow(Nt, d)
            pw = self.fp_pow(self.fpscale(W, fac), 2*k-d)
            term = self.fpmul(term, pw)
            Tt = self.fpadd(Tt, self.fpscale(term, F(tau)))
        Ck = (F(1)+r)**(2*k-1) / (F(2)**(2*k) * (F(1)-r)**(2*k-1))
        num = self.fpscale(Tt, Ck)
        # den = w (w^2+r)^{2k-2} Nt N1t
        N1t = self.fpmul([F(1), F(-1)], [r, F(1)])          # (1-w)(w+r)
        w2r = [r, F(0), F(1)]
        den = list(W)
        den = self.fpmul(den, self.fp_pow(w2r, 2*k-2))
        den = self.fpmul(den, self.fpmul(Nt, N1t))
        return num, den

    def AB_nd(self, k):
        """R^A: (num, den)."""
        F = self.F; r = self.r
        w2r = [r, F(0), F(1)]
        num = self.fpscale(w2r, F(Fr(-1, 2*factorial(k))))
        Nt = self.fpmul([-r, F(1)], [F(1), F(1)])
        den = self.fpmul([F(0), F(1)], Nt)
        return num, den
