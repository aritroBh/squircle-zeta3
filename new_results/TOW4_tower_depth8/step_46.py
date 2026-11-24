work = '<project>/research/08_theory/TOW4_tower_depth8'
script = '''"""prove_adj78.py -- PROOF of the Newton identities at depths 7,8 for ADJACENT keys (p,p-1), all p.
Identical structure to TOW3 prove_adj56.py (proved at delta<=6; byte-inherited machinery):
the off-diagonal e-family is tow4lib.forced (FORCED ascending family, fit-free, PROVED: built from
proved a_delta via the recursion; closing conditions F(delta,floor(delta/2)+1)=0 proved in forced78.py).
tau constants extended to n=8 (tow4lib, verified exactly vs the tow2lib definitions).
Each identity e(delta,s)+e(delta,s+1) == Lambda^{(delta)}(s) is an exact rational identity in p over Q(i)
(diagonal ehat-conversion + fold). Numerator 0 = proof.
Usage: python prove_adj78.py DELTA S
"""
import sys, time
import sympy as sp
import tow4lib as T4

p = sp.Symbol('p')
d = 2*p-1
i = sp.I

def ff(x, c):
    r = sp.Integer(1)
    for u in range(c): r *= (x-u)
    return r

def eR(delta, s, dd, mm, conj=False):
    if s > delta//2: return sp.Integer(0)
    return T4.forced(delta, s, dd, mm, conj)

def tau1(n, conj=False):
    return (-1)**n * T4.tau2f(n) if conj else T4.tau1f(n)
def tau2(n, conj=False):
    return (-1)**n * T4.tau1f(n) if conj else T4.tau2f(n)

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
        terms.append(-p * sp.expand(g) * (-i)**(-4-t) * binom * i**(-6-delta))
    return terms

def main(delta, s):
    lhs = sp.expand(eR(delta, s, d, p) + eR(delta, s+1, d, p))
    terms = feed_terms(delta, s)
    t0 = time.time()
    terms = [sp.cancel(t) for t in terms]
    lhs = sp.cancel(lhs)
    dens = [sp.denom(t) for t in terms] + [sp.denom(lhs)]
    DEN = sp.Integer(1)
    for dn in dens: DEN *= dn
    rem = sp.expand((lhs - sum(terms)) * DEN)
    is0 = sp.simplify(rem) == 0
    print(f'adjacent ({delta},{s}): remainder*DEN == 0 : {is0}  ({round(time.time()-t0,1)}s)', flush=True)
    assert is0

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
'''
open(work + '/prove_adj78.py','w').write(script)
print("wrote prove_adj78.py")
# run (7,3) in foreground with timing (should be fast like (6,3)=0.8s)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run([PY, work+'/prove_adj78.py','7','3'], capture_output=True, text=True, timeout=600)
print(r.stdout, r.stderr[-500:] if r.stderr else '', round(time.time()-t0,1),'s')
