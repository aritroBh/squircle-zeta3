import sys, itertools
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import zsum
from hpl import coeffs_exact

def predict(w):
    nz = [t for t in range(len(w)) if w[t] != 0]
    b = [w[t] for t in nz]; s = len(b)
    m0 = nz[0]
    m = [nz[t+1]-nz[t]-1 for t in range(s-1)]
    eps = b[0]; k = 1+m0; cf = 1
    for x in b: cf *= x
    v = tuple(b[t]*b[t+1]*(1+m[t]) for t in range(s-1))
    return cf, eps, k, v

def inverse(cf, eps, k, v):
    s = len(v)+1; b=[eps]
    for t in range(s-1): b.append(b[-1]*(1 if v[t]>0 else -1))
    w=[0]*(k-1)
    for t in range(s):
        w.append(b[t])
        if t<s-1: w += [0]*(abs(v[t])-1)
    return tuple(w)

# INDEPENDENT ground truth: h_w(N) direct from the HPL integral recursion (hpl.coeffs_exact),
# compared against cf*eps^N/N^k*Z_v(N-1) evaluated by direct nested sums (zsum.Z_val).
ok_map=ok_bij=ok_wt=True; triples=set(); ntest=0
for L in range(1,7):
    for w in itertools.product([0,1,-1], repeat=L):
        if w[-1]==0: continue
        ntest+=1
        cf,eps,k,v = predict(w)
        ce = coeffs_exact(w, 14)   # direct Taylor coefficients of H_w
        for N in range(3,14):
            val = Fr(cf)*Fr(eps**N, N**k)*zsum.Z_val(v, N-1)
            if val != ce[N]:
                ok_map=False; print('TAYLOR MISMATCH', w, N, val, ce[N]); break
        if (eps,k,v) in triples: ok_bij=False; print('COLLISION', w)
        triples.add((eps,k,v))
        cf2=1
        for x in [t for t in w if t!=0]: cf2*=x
        if inverse(cf2,eps,k,v)!=tuple(w): ok_bij=False; print('INV FAIL', w)
        if len(w) != k + sum(abs(t) for t in v): ok_wt=False; print('WEIGHT FAIL', w)
print(f'Lemma B independent check: {ntest} words (len<=6)')
print('  [c^N]H_w == cf eps^N N^-k Z_v(N-1) for N=3..13 (vs DIRECT Taylor):', ok_map)
print('  w->(eps,k,v) injective, inverse recovers w:', ok_bij)
print('  |w| = k+|v| for all:', ok_wt)
# surjectivity: inverse defined on ALL triples since k>=1 and |v_t|>=1
print('  surjectivity: any (eps,k,v) gives m0=k-1>=0, m_t=|v_t|-1>=0 -> valid word. Count check:',
      ntest, '== 2*(3^6-1)/(3-1) =', 2*(3**6-1)//2)
open('progress.log','a').write(f'step02b: Lemma B independent (direct Taylor, {ntest} words len<=6): map={ok_map} bij={ok_bij} |w|=k+|v|={ok_wt}\n')
