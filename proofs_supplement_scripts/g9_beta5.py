# g9_beta5.py  -- Lemma "Beta moments at weight five" (Section g9:sec:order5).
# Formal Beta moments  Ihat(n,k;a,b) in Q[z2,z3,z5]  (z4 := 2 z2^2/5), a+b <= 5.
#  (1) rules: [z5] Ihat = 24 B(n+1,k+1) iff a+b=5, a,b>=1, else 0;  [z2 z3] Ihat = -12 B iff (a,b) in {(3,2),(2,3)}, else 0
#  (2) value check against quadrature of  int t^n (1-t)^k log^a t log^b (1-t) dt  at 40 digits.
import itertools, sympy as sp, mpmath as mp
e, d = sp.symbols('eps delta')
z2, z3, z5 = sp.symbols('z2 z3 z5')
zf = {2: z2, 3: z3, 4: sp.Rational(2, 5)*z2**2, 5: z5}
def H(i, M):            # H_M^{(i)}
    return sum(sp.Rational(1, j**i) for j in range(1, M+1))
def psi(i, M):          # psi^{(i-1)}(M) at an integer M>=1, gamma dropped for i=1
    if i == 1: return H(1, M-1)
    return (-1)**i*sp.factorial(i-1)*(zf[i] - H(i, M-1))
def trunc(p, N=5):
    p = sp.Poly(sp.expand(p), e, d)
    return sum(c*e**i*d**j for (i, j), c in p.terms() if i+j <= N)
def Ihat(n, k, N=5):
    Th = 0
    for i in range(1, N+1):
        Th += sp.Rational(1, sp.factorial(i))*(psi(i, n+1)*e**i + psi(i, k+1)*d**i - psi(i, n+k+2)*(e+d)**i)
    Th = trunc(Th, N)
    ex, pw = sp.Integer(1), sp.Integer(1)
    for j in range(1, N+1):
        pw = trunc(pw*Th, N); ex += pw/sp.factorial(j)
    ex = sp.Poly(sp.expand(ex), e, d)
    Bn = sp.factorial(n)*sp.factorial(k)/sp.factorial(n+k+1)
    return {(a, b): sp.expand(sp.factorial(a)*sp.factorial(b)*Bn*ex.coeff_monomial(e**a*d**b))
            for a in range(N+1) for b in range(N+1-a)}, Bn
def main():
    bad = 0; cnt = 0
    mp.mp.dps = 40
    zv = {z2: mp.zeta(2), z3: mp.zeta(3), z5: mp.zeta(5)}
    worst = mp.mpf(0)
    for n in range(0, 6):
        for k in range(0, 4):
            tab, Bn = Ihat(n, k)
            for (a, b), v in tab.items():
                P = sp.Poly(v, z2, z3, z5)
                c5 = P.coeff_monomial(z5); c23 = P.coeff_monomial(z2*z3)
                w5 = 24*Bn if (a+b == 5 and a >= 1 and b >= 1) else 0
                w23 = -12*Bn if (a, b) in ((3, 2), (2, 3)) else 0
                cnt += 1
                if c5 != w5 or c23 != w23: bad += 1; print('RULE FAILS', n, k, a, b, c5, c23)
                # z5 and z2 z3 must not occur in lower weight; no other monomial contains z5 or z3 z2 beyond these
                for (i, j, l), c in P.terms():
                    if l >= 1 and (i, j, l) != (0, 0, 1) and c != 0: bad += 1; print('z5 times something', n, k, a, b, (i, j, l))
                if (n, k) in ((0, 0), (1, 0), (2, 1), (3, 2), (5, 3)) :
                    val = mp.mpf(0)
                    for (i, j, l), c in P.terms():
                        val += mp.mpf(c.p)/mp.mpf(c.q)*zv[z2]**i*zv[z3]**j*zv[z5]**l
                    num = mp.quad(lambda t: t**n*(1-t)**k*mp.log(t)**a*mp.log(1-t)**b, [0, mp.mpf(1)/2, 1])
                    worst = max(worst, abs(val-num))
    print('entries (n<=5, k<=3, a+b<=5):', cnt, ' rule failures:', bad, ' worst |formal value - quadrature| (5 pairs (n,k), all a+b<=5, 40 digits):', mp.nstr(worst, 3))

if __name__ == '__main__':
    main()
