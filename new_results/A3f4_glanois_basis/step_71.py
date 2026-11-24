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
    print('R%d: residual %s  height 10^%.1f  ratio 10^%.0f' % (
        i+1, mp.nstr(abs(res),3), float(mp.log10(h)), float(mp.log10(abs(res)/h))))

R = sp.Matrix(rels)
snf = R.smith_normal_form()
print('SNF diagonal:', [int(snf[i,i]) for i in range(5)])
