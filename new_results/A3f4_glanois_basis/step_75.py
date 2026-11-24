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

R = sp.Matrix(rels)
print('SNF...')
d = R.smith_normal_form() if hasattr(R, 'smith_normal_form') else None
try:
    snf = sp.matrices.normalforms.smith_normal_form(R)
except Exception as e:
    print('normalforms path:', e)
    snf = None
if snf is not None:
    print('SNF diagonal:', [int(snf[i,i]) for i in range(5)])

# LLL-reduce the 5x32 integer matrix
print('LLL on relation lattice...')
L = R.lll(delta=sp.Rational(99,100))
out = []
for r in L.tolist():
    if all(x==0 for x in r): continue
    res = sum(x*fv for x,fv in zip(r, free_vals))
    h = max(abs(x) for x in r)
    out.append((r,h,res))
    print('height 10^%.1f residual %s' % (float(mp.log10(h)), mp.nstr(abs(res),3)))
print('num vectors:', len(out))
with open(WORK+'/progress.log','a') as f:
    f.write('step_29: SNF + LLL reduction of the 5-dim relation lattice\n')
