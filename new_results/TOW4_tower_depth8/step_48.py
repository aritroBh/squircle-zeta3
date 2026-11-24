work = '<project>/research/08_theory/TOW4_tower_depth8'
script = '''"""prove_adj78_fast.py -- same identities as prove_adj78.py, FAST exact lcm-numerator method.
The identity lhs == sum(terms) is proved by exact rational arithmetic in p over Q(i):
each expression is combined to a single fraction (sp.cancel), denominators are factored into
monic linear factors (all denominators here are products of linear forms (4p-c), p, p+1 and
rational constants), DEN = lcm (product factor^maxexp), and the exact polynomial numerator
    total = lhs_num*(DEN/lhs_den) - sum term_num*(DEN/term_den)
is computed by exact polynomial division and expansion. total == 0 (the zero polynomial) is the proof.
Usage: python prove_adj78_fast.py DELTA S
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

def monic(f):
    f = sp.expand(f)
    if f == 0: return f
    lc = sp.LC(sp.Poly(f, p))
    return sp.expand(f/lc), lc

def main(delta, s):
    t0 = time.time()
    lhs = eR(delta, s, d, p) + eR(delta, s+1, d, p)
    terms = feed_terms(delta, s)
    print(f'  built lhs + {len(terms)} terms ({round(time.time()-t0,1)}s)', flush=True)
    exprs = [lhs] + [-t for t in terms]
    rats = []
    for e in exprs:
        f = sp.cancel(e)
        rats.append(sp.fraction(f))
    print(f'  combined to single fractions ({round(time.time()-t0,1)}s)', flush=True)
    # lcm of denominators via monic linear factor bookkeeping
    fac_pow = {}
    deninfo = []
    for num, den in rats:
        fl = sp.factor_list(den)
        content, flist = fl
        cur = {}
        for f, e in flist:
            fm, lc = monic(f)
            if fm == 1:
                content *= lc
                continue
            fac_pow[fm] = max(fac_pow.get(fm, 0), e)
            cur[fm] = e
        deninfo.append((content, cur))
    DEN = sp.Integer(1)
    for f, e in fac_pow.items():
        DEN *= f**e
    DEN = sp.expand(DEN)
    print(f'  DEN degree {sp.Poly(DEN,p).degree()} from {len(fac_pow)} distinct factors ({round(time.time()-t0,1)}s)', flush=True)
    total = sp.Integer(0)
    for (num, den), (content, cur) in zip(rats, deninfo):
        quot = DEN
        for f, e in cur.items():
            # divide DEN by f^e via exact division
            for _ in range(e):
                q, r = sp.div(quot, f, p)
                assert r == 0, 'nonexact division!'
                quot = q
        total += num * quot / content
    total = sp.expand(total)
    is0 = (total == 0)
    print(f'adjacent ({delta},{s}) FAST: exact numerator == 0 : {is0}  ({round(time.time()-t0,1)}s)', flush=True)
    assert is0

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
'''
open(work + '/prove_adj78_fast.py','w').write(script)
print("wrote prove_adj78_fast.py")
import subprocess, time
PY = '<project>/.venv/bin/python'
for ds in [('7','2'), ('8','3')]:
    t0=time.time()
    r = subprocess.run([PY, work+'/prove_adj78_fast.py']+list(ds), capture_output=True, text=True, timeout=900)
    print(r.stdout.strip(), r.stderr[-400:] if r.stderr else '', 'wall', round(time.time()-t0,1),'s', flush=True)
