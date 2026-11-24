"""prove_E3b56.py -- PROOF of diagonal blindness (E3) at depths delta = 5, 6, all p.
Identical structure to TOW2 prove_E3b34.py (proved at delta 3,4); e-family via tow3lib.
The feed of Delta_pp(delta') = rho_pp(delta') - sum_s (-1)^s ehat_pp(delta',s) (delta' <= delta)
into the child (p,p-1) at depth delta (DL2-direct + fold) vanishes identically.
Usage: python prove_E3b56.py DELTA
"""
import sympy as sp, time, sys
import tow3lib as T3

p = sp.Symbol('p')
i = sp.I

def eR(delta, s, dd, mm, conj=False):
    if s > delta//2: return sp.Integer(0)
    return T3.e(delta, s, dd, mm, conj)

def tau1(n, conj=False):
    return (-1)**n * T3.tau2(n) if conj else T3.tau1(n)
def tau2(n, conj=False):
    return (-1)**n * T3.tau1(n) if conj else T3.tau2(n)

rB0d = 2*i*(4*p-1)/(p+1)

def ehat(delta, s, conj=False):
    tot = 0
    for dpp in range(0, delta+1):
        tot += tau1(delta-dpp, conj) * eR(dpp, s, 2*p+1, p+1, conj) * (rB0d if not conj else -rB0d) / (4*p-dpp-1)
    return sp.expand((p+1)*tot)

def rho_pp(delta, conj=False):
    return T3.e(delta, 0, 2*p, p, conj)

def Delta(delta, conj=False):
    D = rho_pp(delta, conj) - sum( (-1)**s * ehat(delta, s, conj) for s in range(0, delta//2+1) )
    return sp.cancel(D)

def Dc(delta):
    return sp.expand(-Delta(delta, conj=True))

def gfun(t):
    return sum( -Dc(dp)*tau1(t-dp, conj=True)/(4*p-3-dp) for dp in range(0, t+1) )

def feed(delta):
    direct = p*sum( tau2(delta-dp)*Delta(dp)/(4*p-3-dp) for dp in range(0, delta+1) )
    ft = 0
    for t in range(0, delta+1):
        c = (-i)**(-4-t) * i**(-6-delta)
        f = sp.Integer(1)
        for u in range(delta-t): f *= (4*p-6-t-u)
        binom = f/sp.factorial(delta-t)
        ft += -p * sp.expand(gfun(t)) * c * binom
    return sp.expand(direct + ft)

if __name__ == '__main__':
    delta = int(sys.argv[1])
    t0 = time.time()
    f = feed(delta)
    rem = sp.together(f)
    num = sp.expand(sp.fraction(rem)[0])
    is0 = num == 0
    print(f'E3 blindness delta={delta}: numerator == 0 : {is0}  ({round(time.time()-t0,1)}s)', flush=True)
    assert is0
