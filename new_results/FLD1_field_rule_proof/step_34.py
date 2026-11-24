WD='<project>/research/08_theory/FLD1_field_rule_proof'
with open(WD+'/progress.log','a') as fh:
    fh.write("step22: symbolic-n mu_{k,n} (ring in h_r,g_r,z_r over Q(n)), k<=5 matches numeric n<=60; fixed which-bug\n")
code = r'''"""step23_zsum.py -- Z-sum canonical form of mu_{k,n}, per zeta-monomial, general k.

h_r = H^(r)_{n-1} -> 2^{r-1}(Z_r(2n)+Z_{-r}(2n)) - 1/n^r ;  g_r = H^(r)_{2n-1} -> Z_r(2n) - 1/(2n)^r.
Products -> stuffle.  Output Zsum[k][zkey] = {word: coef in Q(n)} (zkey = zeta exponents, padded to 4).
Verified: k=3 reproduces mu_zsum.pkl (from mu3_symbolic, the B3 template) exactly.
"""
import sys, pickle
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0,'.')
import zsum
from zsum import K, nn, ONE, ZERO

src = open("step_22_musymbolic.py").read()
ns = {}
exec(src[:src.index("if __name__")], ns)
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
# rebuild formulas for k<=KMAX (uses ns machinery with its own KMAX=5)
mu_formula = ns['mu_formula']

def harm_n_minus_1(r):
    return {(r,): K(2**(r-1)), (-r,): K(2**(r-1)), (): -ONE/nn**r}
def harm_2n_minus_1(r):
    return {(r,): ONE, (): -ONE/(2*nn)**r}

def to_zsum(F, KMAXf):
    parts = {}
    for mon, cf in F.terms():
        E = {(): ONE}
        for r in range(1, KMAXf+1):
            for _ in range(mon[r-1]):
                E = zsum.mul(E, harm_n_minus_1(r))
            for _ in range(mon[KMAXf+r-1]):
                E = zsum.mul(E, harm_2n_minus_1(r))
        zkey = tuple(list(mon[2*KMAXf:])+[0]*(4-len(mon[2*KMAXf:])))[:4]
        parts.setdefault(zkey, {})
        for w, c in E.items():
            zsum.add_into(parts[zkey], w, cf*c)
    return parts

if __name__ == '__main__':
    Zout = {}
    for k in range(KMAX+1):
        F = mu_formula(k)
        Z = to_zsum(F, ns['KMAX'])
        Zout[k] = Z
        tot_words = sum(len(E) for E in Z.values())
        print(f'k={k}: zeta-monomials={ {zk: len(E) for zk,E in Z.items()} }  total words={tot_words}')
    # verify k=3 against mu_zsum.pkl (B3 template, KMAX=3 there -> zkeys length1,2)
    B='<project>/research/10_frontier/W1_corner_family/B_taylor_pi4/'
    Zmu = pickle.load(open(B+'mu_zsum.pkl','rb'))
    keymap = {'1': (0,0,0,0), 'z2': (1,0,0,0), 'z3': (0,1,0,0)}
    def toKf(s):
        e = sp.sympify(s); num, den = sp.fraction(sp.together(e))
        def polyK(p):
            P = sp.Poly(sp.expand(p), sp.Symbol('n')); out = K(0)
            for (kk,), c in P.terms():
                c = sp.Rational(c); out += K(int(c.p))/K(int(c.q))*nn**kk
            return out
        return polyK(num)/polyK(den)
    ok = True
    for part, zkey in keymap.items():
        muE = {tuple(w): toKf(c) for w, c in Zmu[part]}
        mine = Zout[3].get(zkey, {})
        # compare
        allw = set(muE)|set(mine)
        for w in allw:
            a = muE.get(w, ZERO); b = mine.get(w, ZERO)
            if (a-b) != 0:
                ok=False; print('ZSUM MISMATCH k=3', part, w, 'theirs', a.as_expr() if a else 0, 'mine', b.as_expr() if b else 0)
    print('k=3 Z-sum canonical form reproduces mu_zsum.pkl (B3 template):', ok)
    pickle.dump({k:{zk:[(w, c.as_expr()) for w,c in E.items()] for zk,E in Z.items()} for k,Z in Zout.items()},
                open('mu_zsum_K%d.pkl'%KMAX,'wb'))
    print('saved mu_zsum_K%d.pkl'%KMAX)
'''
open(WD+'/step_23_zsum.py','w').write(code)
import subprocess
VENV='<project>/.venv/bin/python'
r = subprocess.run(['nice',VENV,'step_23_zsum.py','5'],capture_output=True,text=True,cwd=WD)
for line in r.stdout.splitlines():
    if 'ZSUM MISMATCH' in line: print(line[:150]); continue
    print(line)
print(r.stderr[-2000:] if r.returncode else '')