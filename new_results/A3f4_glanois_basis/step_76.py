import os, pickle, time
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
import mpmath as mp
mp.mp.dps = 750
import sympy as sp

with open(WORK+'/vals_hi.pickle','rb') as f: V = pickle.load(f)
free_vals = [mp.mpc(mp.mpf(re), mp.mpf(im)) for re,im in V['vals'][27:]]
freewords = V['freewords']

C = 10**700
m = 32
rows = []
for i,v in enumerate(free_vals):
    r = [0]*(m+2); r[i]=1
    r[m] = int(mp.nint(C*v.real)); r[m+1] = int(mp.nint(C*v.imag))
    rows.append(r)
t=time.time()
R = sp.Matrix(rows).lll(delta=sp.Rational(99,100))
print('lll time', f'{time.time()-t:.0f}s')
short = []
for r in R.tolist():
    tail = max(abs(r[32]), abs(r[33]))
    h = max(abs(x) for x in r[:32])
    short.append((tail, h, r[:32]))
short.sort()
for tail, h, coeffs in short[:8]:
    res = sum(x*fv for x,fv in zip(coeffs, free_vals))
    print('height 10^%.1f tail 10^%.1f resid-ratio 10^%.0f' % (float(mp.log10(h)), float(mp.log10(tail)), float(mp.log10(abs(res)/h))))
with open(WORK+'/progress.log','a') as f:
    f.write('step_30: direct LLL on 32 free values; 5-dim relation lattice shortest vectors\n')
