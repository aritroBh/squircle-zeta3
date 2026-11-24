import sys, os, pickle, time
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 750
import sympy as sp

with open(WORK+'/vals_hi.pickle','rb') as f: V = pickle.load(f)
basis_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][:27]]
free_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][27:]]
freewords = V['freewords']

C = 10**700
def lll_rows(vals):
    m = len(vals)
    rows = []
    for i,v in enumerate(vals):
        r = [0]*(m+2); r[i]=1
        r[m] = int(mp.nint(C*v.real)); r[m+1] = int(mp.nint(C*v.imag))
        rows.append(r)
    R = sp.Matrix(rows).lll(delta=sp.Rational(99,100))
    out = []
    for r in R.tolist():
        out.append(r)
    return out

t=time.time()
R = lll_rows(basis_vals)
tails = []
for r in R:
    tail = max(abs(r[27]), abs(r[28]))
    tails.append(tail)
print('basis-only LLL time', f'{time.time()-t:.0f}s')
print('min tail: 10^%.1f (noise floor C^(2/27) = 10^%.1f)' % (mp.log10(min(tails)), 700*2/27))
# any row with tail << floor would be a relation
for r in R:
    tail = max(abs(r[27]), abs(r[28]))
    if tail < 10**600:
        res = sum(c*v for c,v in zip(r[:27], basis_vals))
        print('CANDIDATE relation, tail 10^%.1f' % mp.log10(tail), 'resid', mp.nstr(abs(res),3))
