"""prove_adj34.py -- PROOF of the Newton identities at depths 3,4 for ADJACENT keys (p,p-1), all p.
Usage: python prove_adj34.py DELTA S
Each identity  e(delta,s)+e(delta,s+1) == Lambda^(delta)(s)  is an exact rational identity in p
over Q(i) (diagonal ehat-conversion + fold included). Remainder 0 = complete proof.
"""
import sys, time
import sympy as sp
import tow2lib as T

p = sp.Symbol('p')
d = 2*p-1
i = sp.I
dd_sym, mm_sym = T.d, T.m

def ff(x, c):
    r = sp.Integer(1)
    for u in range(c): r *= (x-u)
    return r

def eR(delta, s, dd, mm, conj=False):
    if s > delta//2: return sp.Integer(0)
    ii = -sp.I if conj else sp.I
    if s == 0: val = T.arho_simple(delta).subs({dd_sym: dd, mm_sym: mm, sp.I: ii})
    elif delta == 2 and s == 1: val = -1/(2*dd-3)
    elif delta == 3 and s == 1: val = ii*(dd-3)/(2*dd-3)
    elif delta == 4 and s == 1:
        val = (dd**3 - 9*dd**2 + sp.Rational(113,4)*dd - ii*dd/2 + ii*mm - 29)/((2*dd-3)*(2*dd-5))
    elif delta == 4 and s == 2: val = sp.Rational(3)/((2*dd-3)*(2*dd-5))
    else: raise ValueError((delta,s))
    return val

def tau1(n, conj=False):
    return (-1)**n * T.tau2f(n) if conj else T.tau1f(n)
def tau2(n, conj=False):
    return (-1)**n * T.tau1f(n) if conj else T.tau2f(n)

rA0 = i*(4*p-3)/(p+1)
rB0 = (4*p-1)*(4*p-3)/((p+1)*p)

def ehat(delta, s, conj=False):
    tot = 0
    for dpp in range(0, delta+1):
        tot += tau1(delta-dpp, conj) * eR(dpp, s, 2*p+1, p+1, conj) * rB0 / (4*p-dpp-1)
    return sp.expand((p+1)*tot)

def Ediag(delta, s, conj=False):
    tot = 0; t = 0
    while s+t <= delta//2:
        tot += (-1)**t * ehat(delta, s+t, conj)
        t += 1
    return sp.expand(tot)

def feed_terms(delta, s):
    terms = []
    for dp in range(0, delta+1):
        terms.append((p+1)*tau1(delta-dp)*eR(dp, s, 2*p, p+1)*rA0/(2*d-dp-1))
    for dp in range(0, delta+1):
        terms.append(p*tau2(delta-dp)*Ediag(dp, s)/(2*d-dp-1))
    for t in range(0, delta+1):
        mu = 4*p-4-t
        g = 0
        for dp in range(0, t+1):
            g -= Ediag(dp, s, conj=True)*tau1(t-dp, conj=True)/(4*p-3-dp)
        binom = ff(mu-2, delta-t)/sp.factorial(delta-t)
        # (-i)**(4p-4-t) = (-i)**(-4-t) and i**(4p-6-delta) = i**(-6-delta)  [x**4 = 1]
        terms.append(-p * sp.expand(g) * (-i)**(-4-t) * binom * i**(-6-delta))
    return terms

def main(delta, s):
    lhs = sp.expand(eR(delta, s, d, p) + eR(delta, s+1, d, p))
    terms = feed_terms(delta, s)
    t0 = time.time()
    terms = [sp.cancel(t) for t in terms]
    dens = [sp.denom(t) for t in terms] + [sp.denom(lhs)]
    DEN = sp.Integer(1)
    for dn in dens: DEN *= dn
    rem = sp.expand((lhs - sum(terms)) * DEN)
    is0 = sp.simplify(rem) == 0
    print(f'adjacent ({delta},{s}): remainder*DEN == 0 : {is0}  ({round(time.time()-t0,1)}s)', flush=True)
    assert is0

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
