
import sympy as sp

sig = sp.Symbol('sig')
q   = sp.Symbol('q')          # q = sqrt(Z)
ORD = 9
TAU = 12

def T(e, n=ORD):
    return sp.series(e, sig, 0, n).removeO()

y = sum(sp.binomial(sp.Rational(1,2), k)*(q**2*sig**2/4)**k for k in range(0, ORD//2+2)) - q*sig/2
y = T(sp.expand(y))
xs = T(sp.expand(q*y/sig))
Z = q**2
Nsig = sig**-2
def phider(j):
    return T(sp.expand(Nsig*(-1)**(j-1)*sp.factorial(j-1)*((xs+Z)**(-sp.Integer(j)) - xs**(-sp.Integer(j)))))
phi2 = phider(2)
phi_ = {j: phider(j) for j in range(3, TAU+3)}

def lh_coeffs(mval):
    mm = sp.Rational(2*mval+3, 2)
    return [sp.Integer(0)] + [T(sp.expand(sp.Rational((-1)**(j+1), j)*(mm/(xs+Z)**j - 1/xs**j))) for j in range(1, TAU+1)]

def series_exp(coeffs):
    out = [sp.Integer(1)] + [sp.Integer(0)]*TAU
    for n in range(1, TAU+1):
        s = sum(k*coeffs[k]*out[n-k] for k in range(1, n+1))
        out[n] = T(sp.expand(s/n))
    return out

def series_mul(A, B):
    out = [sp.Integer(0)]*(TAU+1)
    for i in range(TAU+1):
        if A[i]==0: continue
        for j in range(TAU+1-i):
            if B[j]==0: continue
            out[i+j] = T(out[i+j] + A[i]*B[j])
    return out

def bracket(mval):
    H = series_exp(lh_coeffs(mval))
    E = [sp.Integer(0)]*(TAU+1)
    for j in range(3, TAU+1):
        E[j] = -phi_[j]/sp.factorial(j)
    ex = series_exp(E)
    integr = series_mul(H, ex)
    br = sp.Integer(0)
    for r in range(0, TAU//2+1):
        br = T(sp.expand(br + sp.expand(integr[2*r]*sp.factorial2(2*r-1)/phi2**r)))
    return sp.expand(br)

# ---- phi_s and Efac ----
logy = T(sp.log(1+(y-1)))
phi_s = T(sp.expand(q*y/sig - 2*logy/sig**2))
Efac = T(sp.exp(T(sp.expand(-phi_s + 2*q/sig - Z/2))))

# ---- Stirling check (numeric) ----
import mpmath
mpmath.mp.dps = 50
for m in [0,1,2]:
    f = lambda N: mpmath.gamma(N-m-mpmath.mpf('1.5'))/mpmath.gamma(N)*N**(m+mpmath.mpf('1.5'))
    c1 = sp.Rational((2*m+3)*(2*m+5),8)
    c2 = sp.Rational((2*m+3)*(2*m+5)*(2*m+7)*(6*m+11),384)
    for N in [2000.0, 5000.0]:
        lhs = f(N); pred = 1 + float(c1)/N + float(c2)/N**2
        print("Stirling m=%d N=%g: exact=%.16f pred=%.16f diff=%.2e" % (m, N, lhs, pred, abs(lhs-pred)))

def series_m(mval):
    ypow = T(sp.expand(y**(-sp.Rational(2*mval+5,2))), 6)
    Gfac = T(sp.expand((y + sig*q/2)**(-sp.Rational(1,2))), 6)
    br = bracket(mval)
    c1 = sp.Rational((2*mval+3)*(2*mval+5),8)
    c2 = sp.Rational((2*mval+3)*(2*mval+5)*(2*mval+7)*(6*mval+11),384)
    S = 1 + c1*sig**2 + c2*sig**4
    tot = T(sp.expand(ypow*Gfac), 6)
    tot = T(sp.expand(tot*br), 6)
    tot = T(sp.expand(tot*Efac), 6)
    tot = T(sp.expand(tot*S), 6)
    return sp.expand(tot)

g1, g2 = sp.symbols('g1 g2')
s0 = series_m(0); s1 = series_m(1); s2 = series_m(2)
print("series_m(0) =", s0)
print("series_m(1) =", s1)
print("series_m(2) =", s2)
full = s0 + g1*sp.Rational(3,2)*q*sig*s1 + g2*sp.Rational(15,4)*q**2*sig**2*s2
full = T(sp.expand(full), 3)
b1 = sp.expand(full.coeff(sig,1))
b2 = sp.expand(full.coeff(sig,2))
print("beta1 =", b1)
print("beta2 =", b2)
