import sys, pickle, json
sys.path.insert(0,'.')
from fractions import Fraction as Fr
from hpl import coeffs_exact

store = pickle.load(open('hpl_constructed_K7.pkl','rb'))
MU = json.load(open('mu_exact_K7_N30.json'))
ok = True
for zkey, terms in store.items():
    zks = ','.join(map(str,zkey))
    words = set(w for (e,w,q) in terms if len(w)>0)
    wc = {w: coeffs_exact(w, 72) for w in words}
    for N in range(-6, 61):
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
print('k=7 exact [c^N]F == mu_{7,n} (all 15 zeta-monomials, N=-6..60, i.e. n<=30; odd/nonpositive==0):', ok)
open('progress.log','a').write(f'step05d: k=7 exact [c^N]F==mu_{{7,n}} for n<=30 (odd & nonpositive parts 0), all 15 zeta-monomials: {ok}\n')
