import os, pickle
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 750
import sympy as sp

with open(WORK+'/vals_hi.pickle','rb') as f: V = pickle.load(f)
free_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][27:]]
freewords = V['freewords']
with open(WORK+'/five_relations.pickle','rb') as f: F5 = pickle.load(f)
rels = F5['relations']

print('numerical verification (values at 800 digits):')
for i, ints in enumerate(rels):
    res = sum(c*fv for c,fv in zip(ints, free_vals))
    h = max(abs(a) for a in ints)
    print(f'R{i+1}: residual {mp.nstr(abs(res),3)}  height 10^{mp.log10(h):.1f}  ratio 10^{mp.log10(abs(res)/h):.0f}')

R = sp.Matrix(rels)
print('SNF...')
snf = R.smith_normal_form()
print('SNF diagonal:', [snf[i,i] for i in range(5)])
