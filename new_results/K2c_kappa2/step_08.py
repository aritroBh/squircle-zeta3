
import sympy as sp

def bracket_for_m(mval, ORD, TAUORD):
    sig = sp.Symbol('sig', positive=True)
    Z   = sp.Symbol('Z')
    tau = sp.Symbol('tau')
    sqZ = sp.sqrt(Z)
    y = 1 - sqZ/2*sig + Z/8*sig**2 - Z**2/128*sig**4 + Z**3/1024*sig**6 - 5*Z**3*sig**8/sp.Rational(32768) + 7*Z**4*sig**10/sp.Rational(262144)
    # extend y properly: y = sqrt(1+Z sig^2/4) - sqZ sig/2 ; recompute series directly
    y = sp.series(sp.sqrt(1+Z*sig**2/4), sig, 0, ORD+1).removeO() - sqZ*sig/2
    xs = sp.series(sp.expand(sqZ*sig**-1*y), sig, 0, ORD).removeO()
    Nsig = sig**-2
    def phider_at_xs(j):
        expr = Nsig*(-1)**(j-1)*sp.factorial(j-1)*((xs+Z)**(-j) - xs**(-j))
        return sp.series(expr, sig, 0, ORD).removeO()
    phi2 = phider_at_xs(2)
    phi_ = {j: phider_at_xs(j) for j in range(3, TAUORD+3)}
    mm = sp.Rational(2*mval+3, 2)
    lh = mm*(sp.log(1 + tau/(xs+Z))) - sp.log(1 + tau/xs)
    lh = sp.series(lh, tau, 0, TAUORD+1).removeO()
    lh = sp.series(lh, sig, 0, ORD).removeO()
    Hfac = sp.series(sp.exp(lh), tau, 0, TAUORD+1).removeO()
    E = sum(-phi_[j]*tau**j/sp.factorial(j) for j in range(3, TAUORD+3))
    expE = sp.series(sp.exp(E), tau, 0, TAUORD+1).removeO()
    integr = sp.series(sp.expand(expE*Hfac), tau, 0, TAUORD+1).removeO()
    integr = sp.expand(integr)
    br = 0
    for pw in range(0, TAUORD+1, 2):
        c = integr.coeff(tau, pw)
        r = pw//2
        br += sp.expand(c*sp.factorial2(2*r-1)/phi2**r)
    br = sp.series(br, sig, 0, 4).removeO()
    return sp.expand(br)

if __name__ == '__main__':
    import sys
    print("=== evolution of the sig^2 coefficient of the m=0 bracket with tau cutoff (ORD=6 as saddle5) ===")
    for TO in [8, 10, 12, 14]:
        br = bracket_for_m(0, 6, TO)
        c2 = sp.expand(br).coeff(sp.Symbol('sig', positive=True), 2)
        print("TAUORD=%2d: b2 = %s" % (TO, c2))
    print()
    print("=== high-order run ORD=10, TAUORD=16: all three brackets ===")
    sig = sp.Symbol('sig', positive=True)
    Z = sp.Symbol('Z')
    for m in [0,1,2]:
        br = bracket_for_m(m, 10, 16)
        print("m=%d bracket = %s" % (m, br))
    print()
    print("truth (numeric, step_07): b2 m=0: -5/8 ; m=1: -7/4+3/(4Z) ; m=2: -27/8+15/(4Z)")
