import os, sys, time, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import mpmath as mp
mp.mp.dps = 250
from fractions import Fraction as F

r2 = mp.sqrt(2); X = mp.log(2); Y = mp.log(1+r2)
b = r2-1

data = pickle.load(open("symbol_fit2.pkl","rb"))
cols = data['cols']
x = pickle.load(open("symbol_fit2_solution.pkl","rb"))['x']
nz = [(label, kind, xv) for (label, kind), xv in zip(cols, x) if xv.p != 0 or xv.q != 0]
print("nonzero solution terms:", len(nz))

z5 = mp.zeta(5); z3 = mp.zeta(3); pi2 = mp.pi**2
def Li_real(n, z): return mp.re(mp.polylog(n, z))
pval = 2+r2
Li5v = {'b': mp.polylog(5,b), '-b': mp.polylog(5,-b), '1/2': mp.polylog(5,mp.mpf('0.5')),
        '1/sqrt2': mp.polylog(5,1/r2), '2-sqrt2': mp.polylog(5,2-r2), 'p': Li_real(5, pval)}
Li4v = {'b': mp.polylog(4,b), '-b': mp.polylog(4,-b), '1/2': mp.polylog(4,mp.mpf('0.5')),
        '1/sqrt2': mp.polylog(4,1/r2), '2-sqrt2': mp.polylog(4,2-r2), 'p': Li_real(4, pval)}
Li3v = {'b': mp.polylog(3,b), '1/2': mp.polylog(3,mp.mpf('0.5')), '1/sqrt2': mp.polylog(3,1/r2)}
Li2v = {'b': mp.polylog(2,b)}
mp.mp.dps = 60
chk = -mp.quad(lambda t: mp.log(abs(1-pval*t))/t, [0, 1/pval, 1], method='tanh-sinh')
print("Re(polylog(2,p)) vs real integral:", mp.nstr(Li_real(2,pval)-chk, 4))
av = 1/pval
nested = mp.quad(lambda t: mp.quad(lambda s: mp.quad(lambda r: mp.log(abs(1-r*av)), [0,s])/(s), [0,t])/(t), [0,1], method='tanh-sinh')
print("Re(polylog(4,p)) vs nested real quad:", mp.nstr(Li_real(4,pval)+nested, 4))
mp.mp.dps = 250

def Li32(x, y, N=1400):
    tot = mp.mpf(0); xp = x
    for n in range(2, N+1):
        if n > 2: xp *= x
        inner = mp.mpf(0); yp = y
        for m in range(1, n):
            if m > 1: yp *= y
            inner += yp/(m*m)
        tot += xp/(n**3)*inner
    return tot
t0=time.time()
L32_b1 = Li32(b, mp.mpf(1)); L32_h1 = Li32(1/r2, mp.mpf(1))
print("Li32 evals in %.1fs" % (time.time()-t0))
pickle.dump({'Li5v':{k:mp.nstr(v,240) for k,v in Li5v.items()},
             'Li4v':{k:mp.nstr(v,240) for k,v in Li4v.items()},
             'Li3v':{k:mp.nstr(v,240) for k,v in Li3v.items()},
             'Li2v':{k:mp.nstr(v,240) for k,v in Li2v.items()},
             'L32_b1': mp.nstr(L32_b1,240), 'L32_h1': mp.nstr(L32_h1,240)}, open("atom_vals.pkl","wb"))
print("saved atom values")
