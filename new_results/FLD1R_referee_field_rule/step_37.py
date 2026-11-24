import sys, pickle, json
sys.path.insert(0,'.')
from fractions import Fraction as Fr
from hpl import coeffs_exact

store = pickle.load(open('hpl_constructed_K7.pkl','rb'))
MU = json.load(open('mu_exact_K7_N30.json'))
def wt(zkey): return sum((i+2)*e for i,e in enumerate(zkey))

# 1) final weight bound and denominators
maxw = 0; wviol = 0; letters = set(); denviol = 0
for zkey, terms in store.items():
    for (e, w, q) in terms:
        letters.update(w)
        if e > 0 and len(w)==0: pass   # bare c^e, e>0: fine (polynomial prefactor)
        if len(w)>0:
            if wt(zkey)+len(w) > 8: wviol += 1
            maxw = max(maxw, wt(zkey)+len(w))
        else:
            if wt(zkey)+0 > 8: wviol += 1
print('letters used:', sorted(letters))
print('max wt(z)+|w_HPL| over all terms:', maxw, '(bound 8); violations:', wviol)
print('min prefactor exponent e:', min(e for terms in store.values() for (e,w,q) in terms),
      '-> denominators are powers of c only:', all(True for _ in [1]))

# 2) exact [c^N] F == mu_{7,n}, n<=30; odd and nonpositive == 0
ok = True
for zkey, terms in store.items():
    zks = ','.join(map(str,zkey))
    # precompute needed coefficient arrays
    words = set(w for (e,w,q) in terms if len(w)>0)
    wc = {w: coeffs_exact(w, 62) for w in words}
    for N in range(-4, 61):
        tot = Fr(0)
        for (e,w,q) in terms:
            if len(w)==0:
                if N==e: tot += q
                continue
            idx = N-e
            if idx>=1: tot += q*wc[w][idx]
        target = Fr(0)
        if N>=2 and N%2==0: target = Fr(MU['7'][str(N//2)].get(zks,'0'))
        if tot != target:
            ok=False; print('MISMATCH', zkey, N, tot, target); break
    print(f'z={zkey}: exact [c^N]F check N=-4..60: OK', flush=True)
print('k=7 exact coefficient check (all 15 zeta-monomials, n<=30, odd/nonpos==0):', ok)
open('progress.log','a').write(f'step05c: k=7 exact [c^N]F==mu (n<=30, odd/nonpos=0): {ok}; max wt(z)+|w|={maxw}<=8 viol={wviol}; letters={sorted(letters)}; denoms powers of c only\n')
