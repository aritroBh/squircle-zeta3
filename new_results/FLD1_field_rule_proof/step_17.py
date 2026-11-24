import os
WD='<project>/research/08_theory/FLD1_field_rule_proof'
code = r'''"""step20_bijection.py -- THE CORE DICTIONARY: HPL word w  <->  (eps, k, v).

CLAIM. For every word w over {0,1,-1} with nonzero last letter, the N-th Taylor coefficient of H_w(c) is
    h_w(N) := [c^N] H_w(c) = cf_w * eps_w^N / N^{k_w} * Z_{v_w}(N-1),
where, writing the nonzero letters of w as b_1..b_s (left to right; b_s = last letter), with
m_0 = #leading zeros and m_t = #zeros between b_t and b_{t+1}:
    eps_w = b_1,   k_w = 1 + m_0,   cf_w = b_1*...*b_s (= +-1),
    v_w = ( b_1 b_2 (1+m_1), ..., b_{s-1} b_s (1+m_{s-1}) ).
The map w -> (eps_w, k_w, v_w) is a BIJECTION from {such w} onto
    {+-1} x Z_{>=1} x {words over Z\{0}}.
We verify this against zsum.hpl_coeff_terms for ALL words of length <= 5, and verify the inverse map.
"""
import sys, itertools
sys.path.insert(0, '.')
from fractions import Fraction as Fr
import zsum

def predict(w):
    nz = [t for t in range(len(w)) if w[t] != 0]
    b = [w[t] for t in nz]
    s = len(b)
    m0 = nz[0]                                  # leading zeros
    m = [nz[t+1]-nz[t]-1 for t in range(s-1)]   # zeros between consecutive nonzero letters
    eps = b[0]; k = 1+m0; cf = 1
    for x in b: cf *= x
    v = tuple(b[t]*b[t+1]*(1+m[t]) for t in range(s-1))
    return cf, eps, k, v

def inverse(cf, eps, k, v):
    # reconstruct w from (eps,k,v); cf is determined, used only as a check
    s = len(v)+1
    b = [eps]
    for t in range(s-1):
        b.append(b[-1]*(1 if v[t]>0 else -1))
    m0 = k-1
    m = [abs(v[t])-1 for t in range(s-1)]
    w = [0]*m0
    for t in range(s):
        w.append(b[t])
        if t < s-1: w += [0]*m[t]
    # check cf
    cf2 = 1
    for x in b: cf2 *= x
    assert cf2 == cf, (cf, cf2)
    return tuple(w)

def words(maxlen):
    for L in range(1, maxlen+1):
        for w in itertools.product([0,1,-1], repeat=L):
            if w[-1] != 0:
                yield w

ok_map = True; ok_inv = True; ntest = 0
triples = set()
for w in words(5):
    ntest += 1
    cf, eps, k, v = predict(w)
    # compare with zsum.hpl_coeff_terms (single term)
    terms = zsum.hpl_coeff_terms(w)
    assert len(terms) == 1, (w, terms)
    (cfz, epsz, kz, vz) = terms[0]
    if not (Fr(cf) == cfz and eps == epsz and k == kz and tuple(v) == tuple(vz)):
        ok_map = False; print('MISMATCH', w, (cf,eps,k,v), (cfz,epsz,kz,vz))
    # bijectivity: triples unique
    if (eps,k,v) in triples:
        ok_map = False; print('COLLISION', w, (eps,k,v))
    triples.add((eps,k,v))
    # inverse
    if inverse(cf, eps, k, v) != tuple(w):
        ok_inv = False; print('INVERSE FAIL', w)
print(f'tested {ntest} words (length<=5): map-from-w matches hpl_coeff_terms and is injective: {ok_map}')
print(f'inverse reconstructs w for all: {ok_inv}')
print()
print('Surjectivity onto {+-1} x Z>=1 x words: for ANY (eps,k,v) the inverse recipe gives a valid word,')
print('since m0=k-1>=0 and m_t=|v_t|-1>=0.  Hence w -> (eps,k,v) is a bijection.  Examples:')
for (eps,k,v) in [(1,1,()),(-1,1,()),(1,2,()),(1,1,(2,)),(-1,1,(3,)),(1,1,(-2,1)),(1,3,(1,)),(-1,2,(-1,))]:
    cf = eps
    for t in range(len(v)): cf *= (1 if v[t]>0 else -1)
    # cf = product of b's = eps * prod sign(v_t); recompute b product:
    b=[eps]
    for t in range(len(v)): b.append(b[-1]*(1 if v[t]>0 else -1))
    cfp=1
    for x in b: cfp*=x
    w = inverse(cfp, eps, k, v)
    print(f'  (eps={eps},k={k},v={v}) -> w={w}, h_w(N) = {cfp}*({eps})^N/N^{k}*Z{list(v)}(N-1)')
'''
open(WD+'/step_20_bijection.py','w').write(code)
import subprocess
VENV='<project>/.venv/bin/python'
r = subprocess.run(['nice',VENV,'step_20_bijection.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else '')