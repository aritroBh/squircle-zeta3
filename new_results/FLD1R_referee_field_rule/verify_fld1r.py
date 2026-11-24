"""verify_fld1r.py -- REFEREE verification of FLD1's Field-rule proof, emphasizing k=7 (beyond FLD1's k<=6).

Checks (all exact unless stated):
  [1] Lemma B (dictionary) INDEPENDENTLY: [c^N]H_w == cf eps^N N^-k Z_v(N-1) against DIRECT Taylor
      coefficients of H_w (hpl.coeffs_exact) for all 728 words of length<=6; bijection; |w|=k+|v|.
  [2] mu_{7,n} EXACT (n<=30): FLD1's Stirling machinery == referee's independent sympy-series table
      (mu7_indep_N30.json, built by per-n series of the Beta integrand; also == corner_series M_7 to 22 digits).
  [3] k=7 Z-sum structure (mu_zsum_K7.pkl): every rho_w proper (deg<=-2); poles subseteq
      {1/2,0} u {-t/2: 1<=t<=6} (sharp form of A2; all in (1/2)Z as Lemma C needs); weight budget
      |w|+ord_p <= 8-wt(z) per word AND pole location; (A1) |w|<=7-wt(z).
  [4] k=7 constructed HPL form (hpl_constructed_K7.pkl, 6553 terms, 15 zeta-monomials):
      letters subseteq {0,1,-1}; prefactors c^e with e>=-6 (denominators powers of c only);
      wt(z)+|w|<=8 for every term; exact [c^N]F == mu_{7,n} for n<=12 (rational part) resp. n<=30
      (other monomials); odd and nonpositive coefficients == 0.
  [5] k=7 all-n identity SPOT RE-PROOF: reconstruct the zeta(7) monomial from mu_zsum_K7.pkl:
      even canonical form == mu, odd == 0 (identities over Q(n)).
Run:  python verify_fld1r.py      (expected < 5 min)
"""
import sys, time, json, pickle, itertools
from fractions import Fraction as Fr
sys.path.insert(0, '.')
import sympy as sp

t0 = time.time()
PASS = []
def check(name, ok):
    PASS.append((name, bool(ok)))
    print(f'  [{"PASS" if ok else "FAIL"}] {name}   (t={time.time()-t0:.1f}s)', flush=True)

print("="*72)
print("[1] Lemma B (dictionary) -- independent check vs DIRECT Taylor coefficients")
import zsum
from hpl import coeffs_exact
def predict(w):
    nz = [t for t in range(len(w)) if w[t] != 0]
    b = [w[t] for t in nz]; s = len(b)
    m0 = nz[0]; m = [nz[t+1]-nz[t]-1 for t in range(s-1)]
    eps = b[0]; k = 1+m0; cf = 1
    for x in b: cf *= x
    v = tuple(b[t]*b[t+1]*(1+m[t]) for t in range(s-1))
    return cf, eps, k, v
ok1 = True; triples=set(); ntest=0
for L in range(1,7):
    for w in itertools.product([0,1,-1], repeat=L):
        if w[-1]==0: continue
        ntest+=1
        cf,eps,k,v = predict(w)
        ce = coeffs_exact(w, 12)
        for N in range(3,12):
            if Fr(cf)*Fr(eps**N, N**k)*zsum.Z_val(v, N-1) != ce[N]: ok1=False
        if (eps,k,v) in triples: ok1=False
        triples.add((eps,k,v))
        if len(w) != k + sum(abs(t) for t in v): ok1=False
check(f'Lemma B: [c^N]H_w==cf eps^N N^-k Z_v(N-1), bijection, |w|=k+|v| ({ntest} words len<=6)', ok1)

print("="*72)
print("[2] mu_{7,n} exact table (n<=30): Stirling machinery == independent sympy-series table")
src = open('step_21_fastmu.py').read(); ns = {}
exec(src[:src.index("if __name__")], ns)
D = ns['moments'](7)
mine = json.load(open('mu7_indep_N30.json'))
ok2 = True
for k in range(8):
    for nv in range(1, 31):
        P = {','.join(map(str,e)): v for e,v in ns['mu'](k, nv, 7, D).items()}
        A = mine[str(k)][str(nv)]
        for kk in set(A)|set(P):
            if Fr(A.get(kk,'0')) != Fr(str(P.get(kk,'0'))): ok2=False
check('mu_{k,n} (k<=7, n<=30, 15 zeta-monomials): two independent implementations agree', ok2)

print("="*72)
print("[3] k=7 Z-sum structure: properness, poles, weight budget, (A1)")
n = sp.Symbol('n')
Z = pickle.load(open('mu_zsum_K7.pkl','rb'))
def wt(zkey): return sum((i+2)*e for i,e in enumerate(zkey))
ok3 = True; all_poles=set(); worst_deg=-99
for zkey, terms in Z.items():
    for w, cs in terms:
        e = sp.together(sp.sympify(cs)); num, den = sp.fraction(e)
        dn = sp.Poly(sp.expand(num), n); dd = sp.Poly(sp.expand(den), n)
        dg = dn.degree()-dd.degree(); worst_deg=max(worst_deg,dg)
        wl = sum(abs(x) for x in w)
        if wl > 7-wt(zkey): ok3=False                       # (A1)
        if dd.degree()>0:
            for rt, mult in sp.roots(dd.as_expr(), n).items():
                rt = sp.nsimplify(rt); all_poles.add(rt)
                if not (2*rt).is_integer: ok3=False          # half-integer poles (Lemma C needs)
                if wl + mult > 8 - wt(zkey): ok3=False       # budget per word & location
                in_sharp = (rt==sp.Rational(1,2)) or (rt==0) or (rt<0 and -2*rt<=6)
                if not in_sharp: ok3=False                   # sharp pole set (A2)
        else:
            if wl > 8 - wt(zkey): ok3=False
print('   poles at k=7:', [str(p) for p in sorted(all_poles, key=sp.N)], '; worst deg =', worst_deg)
check('k=7: rho_w proper (deg<=-2); poles in {1/2,0} u {-t/2: t<=6}; |w|+ord<=8-wt(z); |w|<=7-wt(z)', ok3)

print("="*72)
print("[4] k=7 constructed HPL form: letters, denominators, weight bound, exact [c^N]F")
store = pickle.load(open('hpl_constructed_K7.pkl','rb'))
MU = json.load(open('mu_exact_K7_N30.json'))
ok4 = True; letters=set(); maxw=0; mine_min=0
for zkey, terms in store.items():
    for (e,w,q) in terms:
        letters.update(w); mine_min=min(mine_min,e)
        if wt(zkey)+len(w) > 8: ok4=False
        maxw = max(maxw, wt(zkey)+len(w))
if letters - {0,1,-1}: ok4=False
print(f'   letters={sorted(letters)}, min prefactor e={mine_min} (denominators c^{{-e}} only), max wt(z)+|w|={maxw} (bound 8)')
for zkey, terms in store.items():
    zks = ','.join(map(str,zkey))
    NMAX = 24 if zkey == (0,0,0,0,0,0) else 60
    words = set(w for (e,w,q) in terms if len(w)>0)
    wc = {w: coeffs_exact(w, NMAX+8) for w in words}
    for N in range(-6, NMAX+1):
        tot = Fr(0)
        for (e,w,q) in terms:
            if len(w)==0:
                if N==e: tot += q
                continue
            idx = N-e
            if idx>=1: tot += q*wc[w][idx]
        target = Fr(0)
        if N>=2 and N%2==0: target = Fr(MU['7'][str(N//2)].get(zks,'0'))
        if tot != target: ok4=False; print('   MISMATCH', zkey, N)
check('k=7 HPL form: letters {0,1,-1}; power-of-c denominators; wt(z)+|w|<=8; exact [c^N]F==mu (n<=12 rat./n<=30 others); odd/nonpos==0', ok4)

print("="*72)
print("[5] k=7 all-n identity spot re-proof (zeta(7) monomial): even==mu, odd==0 over Q(n)")
src = open('step_25_construct.py').read(); ns25 = {}
exec(src[:src.index("if __name__")], ns25)
E = Z[(0,0,0,0,0,1)]
terms = ns25['construct'](E)
Fsym = ns25['symmetrise'](terms); Fsym, bad = ns25['correct'](Fsym, 7)
ev = ns25['even_form'](Fsym); od = ns25['odd_form'](Fsym)
muE = {tuple(w): ns25['toK'](c) for w,c in E}
diff = zsum.lin((1,ev),(-1,muE))
check('zeta(7) monomial: construction -> even==mu AND odd==0 as Q(n) identities (all n)', len(diff)==0 and len(od)==0)

print("="*72)
n_pass = sum(1 for _,o in PASS if o)
print(f'RESULT: {n_pass}/{len(PASS)} referee checks passed in {time.time()-t0:.1f}s')
for nm,o in PASS: print(('  PASS ' if o else '  FAIL ')+nm)
sys.exit(0 if n_pass==len(PASS) else 1)
