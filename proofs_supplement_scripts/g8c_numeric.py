"""g8c_numeric.py -- numerical checks (NOT a proof) for the section 'The weight-three part as an iterated integral'.
At each angle alpha (c=cos(alpha/2), s=sin(alpha/2), k=cos(alpha), Y=1/s), 40 digits:
 (a) M2def   = M_2 from the DEFINITION: Cauchy integral in u (radius 0.2, 64 nodes) of the integral M(u;alpha).
 (b) M2g     = int_0^1 G_2(g) dg, G_2 = (1/g)[RA L^2 + RB L m + RC m^2]   (Lemma 1 and the hand formula).
 (c) M2t     = int_1^{t1} [RA^t L^2 + RB^t L m + RC^t m^2] dt, RA^t = RA(g(t)) dg/dt etc. (chain rule, numerically).
 (d) W2int   = -int_1^{t1} log g dt/(t+1)   versus   zeta(2)/2 - Li2(c^2)/2.
 (e) W3y     = sum_j int_0^Y [A_j L^2 + B_j L m + C_j m^2] dy/(y-a_j)   (the theorem), with L, m from the lambda_j.
 (f) W3     = M2def - W2   versus  W3y, and versus  Y(c)+Y(-c) - zeta(2)(1-c^2)/(2c) artanh c  (HPL series; the
     closed form of Theorem 20), and the sum W2 + W3y versus M2def.
"""
import mpmath as mp
import sys
mp.mp.dps = 48

def M2_def(alpha, r=mp.mpf('0.2'), Nn=64):
    ca = mp.cos(alpha)
    def Mu(u):
        def f(g):
            h = g**(1 - u)
            return ((1 + g)**(-1 - u)*mp.sqrt(1 + g**(2 - 2*u) - 2*ca*h) - 1)*g**(u - 1)
        return mp.quad(f, [0, mp.mpf(1)/2, 1])
    tot = 0
    for j in range(Nn):
        u = r*mp.expjpi(2*mp.mpf(j)/Nn)
        tot += Mu(u)*u**(-2)
    return mp.re(tot/Nn)

def parts(alpha):
    c, s, k = mp.cos(alpha/2), mp.sin(alpha/2), mp.cos(alpha)
    Y = 1/s
    # (b) g-integral
    def G2(g):
        L, m = mp.log(g), mp.log1p(g)
        S = 1 + g*g - 2*k*g
        R = mp.sqrt(S); Rp = (g - k)/R; Rpp = (1 - k*k)/R**3
        RA = ((R - g*Rp + g*g*Rpp)/(2*(1 + g)) - mp.mpf(1)/2)/g
        RB = (g*Rp - R)/(g*(1 + g))
        RC = R/(2*g*(1 + g))
        return RA*L*L + RB*L*m + RC*m*m
    M2g = mp.quad(G2, [0, mp.mpf(1)/16, mp.mpf(1)/2, 1])
    # (c) t-integral (chain rule)
    t1 = 1 + 2*s
    def gt(t): return (t*t - 1)/(2*(t - k))
    def Ft(t):
        g = gt(t); dg = (t*t - 2*k*t + 1)/(2*(t - k)**2)
        return G2(g)*g*dg/g if False else G2(g)*dg
    M2t = mp.quad(Ft, [1, 1 + (t1 - 1)/64, 1 + (t1 - 1)/4, t1])
    # (d)
    W2int = -mp.quad(lambda t: mp.log(gt(t))/(t + 1), [1, 1 + (t1 - 1)/64, 1 + (t1 - 1)/4, t1])
    W2 = (mp.zeta(2) - mp.polylog(2, c*c))/2
    # (e) y-integral
    a = [mp.mpf(0), -1, -1/s**2, -1/(1 + c), -1/(1 - c)]
    qq = (c*c + 1)/(4*c)
    tab = [(0, -1, mp.mpf(1)/2), (mp.mpf(1)/2, 0, mp.mpf(1)/2), (-1, 1, -mp.mpf(1)/2), (-qq, c, -c), (qq, -c, c)]
    def f(y):
        l0, l1, l2, l3, l4 = mp.log(y), mp.log1p(y), mp.log1p(s*s*y), mp.log1p((1 + c)*y), mp.log1p((1 - c)*y)
        L, m = l0 + l2 - l1, l3 + l4 - l1
        return sum((A*L*L + B*L*m + C*m*m)/(y - aj) for (A, B, C), aj in zip(tab, a))
    W3y = mp.quad(f, [0, mp.mpf(1)/1024, mp.mpf(1)/32, mp.mpf(1)/4, 1, Y])
    return M2g, M2t, W2int, W2, W3y

# HPL values from Taylor series (own code, same recursion as Definition of H_w)
def hpl_coeffs(w, N):
    if len(w) == 1:
        a = w[0]
        return [mp.mpf(0)] + [mp.mpf(1)/m if a == 1 else mp.mpf((-1)**(m + 1))/m for m in range(1, N + 1)]
    inner = hpl_coeffs(w[1:], N)
    a = w[0]
    out = [mp.mpf(0)]*(N + 1)
    if a == 0:
        for m in range(1, N + 1): out[m] = inner[m]/m
    else:
        acc = mp.mpf(0); part = []
        for m in range(N + 1):
            acc = acc*a + inner[m]; part.append(acc)
        for m in range(1, N + 1): out[m] = part[m - 1]/m
    return out

def Hval(w, x, N):
    co = hpl_coeffs(w, N)
    return sum(co[m]*x**m for m in range(1, N + 1))

def nz(w): return sum(1 for a in w if a)

def Ysum(c, N):
    def H(w, x):
        return Hval(w, x, N) if x > 0 else (-1)**nz(w)*Hval(tuple(-a for a in w), -x, N)
    def Yf(x):
        return ((1 - x)*(H((0, 0, 1), x) + H((0, 1, -1), x) - H((1, -1, 1), x))
                + (1 - x**2)/(2*x)*H((1, 0, 1), x) - (1 - x)**2/(2*x)*H((1, 0, -1), x))
    return Yf(c) + Yf(-c)

if __name__ == '__main__':
    angles = [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3), ('1.8', mp.mpf('1.8')), ('5pi/6', 5*mp.pi/6)]
    for name, al in angles:
        c = mp.cos(al/2)
        N = int(mp.ceil(50/(-mp.log10(c)))) + 5
        M2d = M2_def(al)
        M2g, M2t, W2int, W2, W3y = parts(al)
        W3 = M2d - W2
        W3cl = Ysum(c, N) - mp.zeta(2)*(1 - c*c)/(2*c)*mp.atanh(c)
        print(f'alpha = {name}  c = {mp.nstr(c, 12)}  s = {mp.nstr(mp.sin(al/2), 12)}  (HPL series order {N})')
        print(f'   M_2 (definition)       = {mp.nstr(M2d, 30)}')
        print(f'   |M2def - M2g|          = {mp.nstr(abs(M2d - M2g), 3)}     (integral of the u^2 coefficient, Lemma 1)')
        print(f'   |M2g   - M2t|          = {mp.nstr(abs(M2g - M2t), 3)}     (Euler substitution)')
        print(f'   |W2int - W2|           = {mp.nstr(abs(W2int - W2), 3)}     (W_2 = zeta(2)/2 - Li2(c^2)/2)')
        print(f'   W_3 = M2def - W2       = {mp.nstr(W3, 30)}')
        print(f'   |W_3 - W3y|            = {mp.nstr(abs(W3 - W3y), 3)}     (the theorem: iterated-integral form)')
        print(f'   |W_3 - (Y(c)+Y(-c)-zeta(2)(1-c^2)/(2c) artanh c)| = {mp.nstr(abs(W3 - W3cl), 3)}')
        print(f'   |W2 + W3y - M2def|     = {mp.nstr(abs(W2 + W3y - M2d), 3)}')
        sys.stdout.flush()
