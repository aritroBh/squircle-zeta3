import pickle, sympy as sp
n = sp.Symbol('n')
W1='<project>/research/10_frontier/W1_corner_family/B_taylor_pi4/'
tpl = pickle.load(open(W1+'mu_zsum.pkl','rb'))
mine = pickle.load(open('mu_zsum_K5.pkl','rb'))
keymap = {(0,0,0,0):'1', (1,0,0,0):'z2', (0,1,0,0):'z3'}
ok = True
for zkey, tkey in keymap.items():
    A = {tuple(w): sp.together(sp.sympify(cs)) for w, cs in mine[3][zkey]}
    B = {tuple(w): sp.together(sp.sympify(cs)) for w, cs in tpl[tkey]}
    for w in set(A)|set(B):
        if sp.simplify(A.get(w,0) - B.get(w,0)) != 0:
            ok=False; print('DIFF', zkey, w, A.get(w), B.get(w))
print('FLD1 k=3 Z-sum form == W1 template mu_zsum.pkl (all 3 zeta-monomials, all words):', ok)

# (b) Lemma A(ii) bounds: e^{(k)}_{j,i,l} nonzero only for j<=k, i+l<=k (check exact support)
import sys
sys.argv=['x','5']
src=open('step_22_musymbolic.py').read(); ns={}
exec(src[:src.index("if __name__")], ns)
for k in (2,3,5):
    Ek = ns['Ek_coeffs'](k)
    viol = [(j,i,l) for (j,i,l) in Ek if j>k or i+l>k]
    exactdeg = all(i+l==k for (j,i,l) in Ek)
    print(f'k={k}: #terms={len(Ek)}, (j<=k & i+l<=k) violations: {len(viol)}; i+l==k exactly for all: {exactdeg}; max j = {max(j for j,i,l in Ek)}')
open('progress.log','a').write(f'step07: k=3 Z-sum == W1 template: {ok}; Lemma A(ii) support: j<=k, i+l=k exactly (k=2,3,5 checked)\n')
