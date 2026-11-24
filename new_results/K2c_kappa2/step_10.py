
import sympy as sp

sig = sp.Symbol('sig')
q   = sp.Symbol('q')          # q = sqrt(Z)
ORD = 9
TAU = 12

def T(e):
    return sp.series(e, sig, 0, ORD).removeO()

# exact-rational y series
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
    return [sp.Integer(0)] + [T(sp.expand((-1)**(j+1)/j*(mm/(xs+Z)**j - 1/xs**j))) for j in range(1, TAU+1)]

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

def integrand_coeffs(mval):
    H = series_exp(lh_coeffs(mval))
    E = [sp.Integer(0)]*(TAU+1)
    for j in range(3, TAU+1):
        E[j] = -phi_[j]/sp.factorial(j)
    ex = series_exp(E)
    return series_mul(H, ex)

def bracket(mval, taumax):
    integr = integrand_coeffs(mval)
    br = sp.Integer(0)
    for r in range(0, taumax//2+1):
        br = T(sp.expand(br + sp.expand(integr[2*r]*sp.factorial2(2*r-1)/phi2**r)))
    return sp.expand(br)

if __name__ == '__main__':
    print("=== m=0 bracket sig^2 coefficient vs tau cutoff (exact) ===")
    for tmax in [8, 10, 12]:
        br = bracket(0, tmax)
        c2 = sp.expand(br).coeff(sig, 2)
        print("taumax=%2d: b2 = %s" % (tmax, sp.simplify(c2)))
    print()
    for m in [0,1,2]:
        br = bracket(m, 12)
        c2 = sp.expand(br).coeff(sig, 2)
        c1 = sp.expand(br).coeff(sig, 1)
        print("m=%d  b1 = %s ; b2 = %s" % (m, sp.simplify(c1), sp.simplify(c2)))
    print()
    print("=== missing pieces: r=5 and r=6 contributions to sig^2 coefficient ===")
    for m in [0,1,2]:
        integr = integrand_coeffs(m)
        for r in [5,6]:
            term = sp.expand(integr[2*r]*sp.factorial2(2*r-1)/phi2**r)
            term = sp.expand(sp.series(term, sig, 0, 3).removeO())
            c2 = term.coeff(sig, 2)
            print("m=%d r=%d: sig^2 part = %s" % (m, r, sp.simplify(c2)))
