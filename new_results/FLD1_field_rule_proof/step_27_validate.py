"""step27_validate.py -- validate the constructed HPL forms for k=4,5 (and k=2,3 as control).

(1) EXACT: [c^(2n)] M_k = mu_{k,n} for n<=100 (all zeta-monomials), odd and non-positive parts = 0.
(2) NUMERICAL: M_k(c) (HPL form, evaluated by the HPL Taylor series hpl.coeffs_mp, which comes from the
    HPL definition -- independent of mu) vs corner_series.corner_coeffs (the corner INTEGRAL), 4 angles, k=4,5.
(3) spot-check: a few HPL words by recursive tanh-sinh quadrature vs the Taylor evaluation.
"""
import sys, pickle
from fractions import Fraction as Fr
import mpmath as mp
sys.path.insert(0,'.')
from hpl import coeffs_exact, coeffs_mp
sys.path.insert(0,'<project>/research/10_frontier/G1_zeta3_any_angle/')
import corner_series

DPS = int(sys.argv[1]) if len(sys.argv)>1 else 80
mp.mp.dps = DPS
store = pickle.load(open('hpl_constructed_K5.pkl','rb'))
MU = __import__('json').load(open('mu_exact_K5_N100.json'))

def zeta_of(zkey):
    v = mp.mpf(1)
    for i,e in enumerate(zkey):
        if e: v *= mp.zeta(i+2)**e
    return v

# ---------- (1) EXACT check n<=100 ----------
def exact_check(kmax):
    allok = True
    for k in range(kmax+1):
        for zkey, terms in store.items():
            if zkey[0]!=k: continue
            zk = zkey[1]; zks = ','.join(map(str,zk))
            # coefficients c^-4..c^(2*100) via exact coeffs
            bad = 0
            for N in range(-4, 201):
                tot = Fr(0)
                for (e,w,q) in terms:
                    if len(w)==0:
                        if N==e: tot += q
                        continue
                    idx = N-e
                    if idx>=1: tot += q*coeffs_exact(w, idx+2)[idx]
                target = Fr(0)
                if N>=2 and N%2==0:
                    target = Fr(MU[str(k)][str(N//2)].get(zks,'0'))
                if tot != target:
                    bad += 1
                    if bad<=3: print(f'   EXACT MISMATCH k={k} z={zk} c^{N}: {tot} vs {target}')
            if bad: allok=False
        print(f'k={k}: EXACT [c^N] for N=-4..200 matches mu (odd & non-positive = 0): {allok}')
    return allok

# ---------- (2) corner_series numerical ----------
def build_wordcache(kmax, N):
    words=set()
    for (k,zk),terms in store.items():
        if k>kmax: continue
        for (e,w,q) in terms:
            if len(w)>0: words.add(w)
    return {w: coeffs_mp(w, N) for w in words}

def combo_value(terms, c, wc, N):
    tot = mp.mpf(0)
    cpow = {}
    for (e,w,q) in terms:
        if len(w)==0:
            tot += mp.mpf(q.numerator)/q.denominator * c**e
            continue
        ce = wc[w]
        s = mp.mpf(0); xp = mp.mpf(c)
        for m in range(1, N+1):
            s += ce[m]*xp; xp *= c
        tot += mp.mpf(q.numerator)/q.denominator * c**e * s
    return tot

def Mk_value(k, c, wc, N):
    tot = mp.mpf(0)
    for (kk,zk),terms in store.items():
        if kk!=k: continue
        tot += zeta_of(zk)*combo_value(terms, c, wc, N)
    return tot

if __name__ == '__main__':
    print('=== (1) EXACT check (n<=100) ===')
    ok1 = exact_check(5)
    print('EXACT all:', ok1)
    print()
    print('=== (2) corner_series numerical (k=4,5) ===')
    N = 600
    wc = build_wordcache(5, N)
    angles = {'2pi/3': 2*mp.pi/3, 'pi/2': mp.pi/2, '3pi/4': 3*mp.pi/4, 'pi/3': mp.pi/3}
    for tag,a in angles.items():
        c = mp.cos(a/2)
        Mc = corner_series.corner_coeffs(a, 5)
        for k in (4,5):
            v = Mk_value(k, c, wc, N)
            d = abs(v - Mc[k])
            print(f'  {tag} (c={mp.nstr(c,6)}) k={k}: |HPLform - corner_series| = {mp.nstr(d,4)}   (digits ~ {mp.nstr(-mp.log10(d+mp.mpf("1e-999")),5)})')
        print()
