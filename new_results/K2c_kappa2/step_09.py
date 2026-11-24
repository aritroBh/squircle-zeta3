
import sympy as sp

sig = sp.Symbol('sig')
q   = sp.Symbol('q')          # q = sqrt(Z)
ORD = 9                        # keep sig-powers < ORD
TAU = 12                       # tau moments up to 2r = TAU

def T(e):
    """truncate sig-series below ORD"""
    return sp.series(e, sig, 0, ORD).removeO()

# ---- saddle x_s : y = sqrt(1+q^2 sig^2/4) - q sig/2 ; x_s = q y / sig ; Z = q^2 ----
y = T(sp.sqrt(1 + q**2*sig**2/4)) - q*sig/2
y = T(sp.expand(y))
xs = T(sp.expand(q*y/sig))                       # x_s
Z = q**2

Nsig = sig**-2
def phider(j):
    return T(sp.expand(Nsig*(-1)**(j-1)*sp.factorial(j-1)*((xs+Z)**(-sp.Integer(j)) - xs**(-sp.Integer(j)))))
phi2 = phider(2)
phi_ = {j: phider(j) for j in range(3, TAU+3)}

# ---- Hfac: exp(lh), lh = mm log(1+tau/(xs+Z)) - log(1+tau/xs), as tau-polynomial with sig-series coeffs ----
def lh_coeffs(mval):
    mm = sp.Rational(2*mval+3, 2)
    coeffs = [sp.Integer(0)]  # tau^0
    for j in range(1, TAU+1):
        c = T(sp.expand((-1)**(j+1)/j * (mm/(xs+Z)**j - 1/xs**j)))
        coeffs.append(c)
    return coeffs

def series_exp(coeffs):
    """given L(tau)=sum c_j tau^j (c_0=0), return coeffs of exp(L)"""
    out = [sp.Integer(1)] + [sp.Integer(0)]*TAU
    for n in range(1, TAU+1):
        s = 0
        for k in range(1, n+1):
            s = s + k*coeffs[k]*out[n-k]
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
    lh = lh_coeffs(mval)
    H = series_exp(lh)
    E = [sp.Integer(0)]*(TAU+1)
    for j in range(3, TAU+1):
        E[j] = -phi_[j]/sp.factorial(j)
    ex = series_exp(E)
    integr = series_mul(H, ex)
    br = sp.Integer(0)
    for r in range(0, TAU//2+1):
        c = integr[2*r]
        br = br + sp.expand(c*sp.factorial2(2*r-1)/phi2**r)
        br = T(sp.expand(br))
    return sp.expand(br)

if __name__ == '__main__':
    import time
    for m in [0,1,2]:
        t0=time.time()
        br = bracket(m)
        print("m=%d bracket = %s   [%.1fs]" % (m, br, time.time()-t0))
    print("truth: b2 m=0: -5/8 ; m=1: -7/4+3/(4Z) ; m=2: -27/8+15/(4Z)  (Z=q^2)")
