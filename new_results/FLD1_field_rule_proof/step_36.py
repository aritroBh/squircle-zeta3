WD='<project>/research/08_theory/FLD1_field_rule_proof'
with open(WD+'/progress.log','a') as fh:
    fh.write("step24: rho_w(n) PROPER (deg-2), poles in {-(k-1)/2..0,1/2} half-integers; no c d/dc needed; only powers of c as denominators\n")
code = r'''"""step25_construct.py -- CONSTRUCT the HPL form from the Z-sum canonical form (the general algorithm).

ALGORITHM (proves Z-sum expression  ==>  HPL combination, even case).  Input: E = {word: rho(n)} with rho in
Q(n) proper, half-integer poles, meaning mu = sum_w rho(n) Z_w(2n).  Output: even HPL combination F with
[c^(2n)]F = E, [c^(2n+1)]F = 0.
  recursion on max Z-sum depth: for the leading word v with coefficient rho(n), partial-fraction rho in m=2n,
  rho = sum_{e,k} c_{e,k}/(2n-e)^k ; each (e,k) gets the HPL word w0 = BJI(+1,k,v) (bijection inverse) and term
  q c^e H_{w0}, q = c_{e,k}/cf.  Since [c^(2n)] c^e H_{w0} = cf/(2n-e)^k Z_v(2n) + (strict suffixes, lower depth),
  the Z_v terms cancel exactly; recurse on the remainder (strictly lower depth).  Base v=(): single-letter HPLs.
  Then SYMMETRISE: F_sym = (F(c)+F(-c))/2 via H_w(-c) = (-1)^{#nonzero letters} H_{nu(w)}(c).
Verification (the PROOF): even canonical Z-sum form of F_sym == E identically in n; odd form == 0; and
low Taylor orders c^-2..c^3 match the exact mu table (step_21).
"""
import sys, pickle
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0,'.')
import zsum
from zsum import K, nn, ONE, ZERO, add_into, lin
from hpl import coeffs_exact

n = sp.Symbol('n')
Zin = pickle.load(open('mu_zsum_K5.pkl','rb'))
MU = __import__('json').load(open('mu_exact_K5_N100.json'))

def toK(s):
    e = sp.together(sp.sympify(s)); num, den = sp.fraction(e)
    def polyK(p):
        P = sp.Poly(sp.expand(p), n); out = K(0)
        for (kk,), c in P.terms():
            c = sp.Rational(c); out += K(int(c.p))/K(int(c.q))*nn**kk
        return out
    return polyK(num)/polyK(den)

def sgn(a): return 1 if a>0 else -1

def bji(eps, k, v):
    """inverse bijection: (eps,k,v) -> HPL word w0 (and cf = product of nonzero letters)."""
    s = len(v)+1
    b = [eps]
    for t in range(s-1): b.append(b[-1]*sgn(v[t]))
    m0 = k-1; m = [abs(v[t])-1 for t in range(s-1)]
    w0 = [0]*m0
    for t in range(s):
        w0.append(b[t])
        if t < s-1: w0 += [0]*m[t]
    cf = 1
    for x in b: cf *= x
    return tuple(w0), cf

def partfrac(sigK):
    """rho(n) in K -> list of (e, k, c_{e,k}) with rho = sum c_{e,k}/(2n-e)^k.  (proper, integer poles in m=2n)"""
    sig = sp.together(sp.sympify(sigK))
    m = sp.Symbol('m')
    R = sp.expand(sig.subs(n, m/2))
    R = sp.together(R)
    out = []
    ap = sp.apart(R, m)
    for term in sp.Add.make_args(ap):
        # each term c/(m-e)^k  (or polynomial, should not occur)
        num, den = sp.fraction(term)
        dd = sp.Poly(sp.expand(den), m)
        if dd.degree() == 0:
            # polynomial part -> should not happen (proper); flag
            out.append(('poly', 0, sp.Rational(num)/sp.Rational(den)))
            continue
        # factor den as (m-e)^k
        roots = sp.roots(sp.expand(den), m)
        assert len(roots) == 1, ('nonlinear den', term)
        e = list(roots.keys())[0]; k = roots[e]
        assert e.is_integer, ('noninteger pole', e)
        c = sp.Rational(num)/(sp.expand(den)/(m-e)**k)
        out.append((int(e), int(k), sp.Rational(c)))
    return out

def construct(E):
    """E = {word: K}  ->  list of (e, w, q Fraction) with even part E (before symmetrisation)."""
    R = {tuple(w): toK(c) for w, c in E}
    terms = []
    guard = 0
    while R and guard < 10000:
        guard += 1
        # pick a max-depth word
        w = max(R.keys(), key=lambda x: (len(x), x))
        sigK = R[w]
        pf = partfrac(sigK.as_expr())
        contrib = {}
        for (e, k, c) in pf:
            if e == 'poly':
                raise ValueError('polynomial part in rho!')
            w0, cf = bji(1, k, w)
            q = Fr(c.numerator, c.denominator)/cf
            terms.append((e, w0, q))
            # canonical form of q * c^e H_{w0} at even part: q * hpl_coeff_at(w0, -e)
            X = zsum.hpl_coeff_at(w0, -e)
            qK = K(q.numerator)/K(q.denominator)
            for ww, cc in X.items():
                add_into(contrib, ww, qK*cc)
        # subtract
        R = lin((1, R), (-1, contrib))
    assert guard < 10000
    return terms

def nu(w):
    return tuple(0 if x==0 else -x for x in w)

def symmetrise(terms):
    """(F(c)+F(-c))/2 as HPL combination over c."""
    out = {}
    for (e, w, q) in terms:
        # F term q c^e H_w(c)
        key=(e,w); out[key]=out.get(key,Fr(0))+q/2
        # F(-c): q (-c)^e H_w(-c) = q (-1)^e (-1)^{s_w} c^e H_{nu(w)}(c)
        s = sum(1 for x in w if x!=0)
        sgnfac = (-1)**e * (-1)**s
        key=(e,nu(w)); out[key]=out.get(key,Fr(0))+q*sgnfac/2
    return [(e,w,q) for (e,w),q in out.items() if q!=0]

def even_form(terms):
    E = {}
    for (e, w, q) in terms:
        X = zsum.hpl_coeff_at(w, -e)
        qK = K(q.numerator)/K(q.denominator)
        for ww, cc in X.items():
            add_into(E, ww, qK*cc)
    return E
def odd_form(terms):
    E = {}
    for (e, w, q) in terms:
        X = zsum.hpl_coeff_at(w, 1-e)
        qK = K(q.numerator)/K(q.denominator)
        for ww, cc in X.items():
            add_into(E, ww, qK*cc)
    return E

def low_order(terms, k, zkey):
    zkeystr = ','.join(map(str,zkey))
    ok=True
    for N in range(-2, 4):
        tot = Fr(0)
        for (e, w, q) in terms:
            idx = N - e
            if idx >= 1:
                tot += q*coeffs_exact(w, idx+2)[idx]
        target = Fr(0)
        if N>=2 and N%2==0:
            target = Fr(MU[str(k)][str(N//2)].get(zkeystr,'0'))
        if tot != target:
            ok=False; print(f'   LOW-ORDER MISMATCH k={k} {zkey} c^{N}: {tot} vs {target}')
    return ok

if __name__ == '__main__':
    KMAX = int(sys.argv[1]) if len(sys.argv)>1 else 3
    store = {}
    for k in range(0, KMAX+1):
        for zkey, E in sorted(Zin[k].items()):
            terms = construct(E)
            Fsym = symmetrise(terms)
            ev = even_form(Fsym)
            od = odd_form(Fsym)
            # even == E ?
            muE = {tuple(w): toK(c) for w,c in E}
            diff = lin((1,ev),(-1,muE))
            lo = low_order(Fsym, k, zkey)
            print(f'k={k} z={zkey}: constructed {len(terms)} terms -> sym {len(Fsym)}; even==mu: {len(diff)==0}; odd==0: {len(od)==0}; low-orders: {lo}')
            if diff: print('   EVEN residual:', list(diff.items())[:3])
            if od: print('   ODD residual:', list(od.items())[:3])
            store[(k,zkey)] = Fsym
    pickle.dump({(k,zk): v for (k,zk),v in store.items()}, open('hpl_constructed_K%d.pkl'%KMAX,'wb'))
    print('saved hpl_constructed_K%d.pkl'%KMAX)
'''
open(WD+'/step_25_construct.py','w').write(code)
import subprocess
VENV='<project>/.venv/bin/python'
r = subprocess.run(['nice',VENV,'step_25_construct.py','3'],capture_output=True,text=True,cwd=WD)
for line in r.stdout.splitlines():
    if 'residual' in line or 'MISMATCH' in line: print(line[:200]); continue
    print(line)
print(r.stderr[-3000:] if r.returncode else '')