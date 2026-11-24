import os, sys, time, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 250
from fractions import Fraction as F

r2 = mp.sqrt(2); X = mp.log(2); Y = mp.log(1+r2)
b = r2-1
def q2f(pq): return mp.mpf(str(pq.p)) + mp.mpf(str(pq.q))*r2

# solution coefficients
data = pickle.load(open("symbol_fit2.pkl","rb"))
cols = data['cols']
x = pickle.load(open("symbol_fit2_solution.pkl","rb"))['x']
terms = [(label, kind, xv) for (label, kind), xv in zip(cols, x) if xv != mp.mpf(0) and str(xv) != '0']
nz = [(label, kind, xv) for (label, kind), xv in zip(cols, x) if xv.p != 0 or xv.q != 0]
print("nonzero solution terms:", len(nz))

# evaluate atom values
z5 = mp.zeta(5); z3 = mp.zeta(3); pi2 = mp.pi**2
def Lin(n, z): return mp.polylog(n, z)
pval = 2+r2
# real-branch Li_n(p) = Re(polylog) -- cross-check below
def Li_real(n, z):
    return mp.re(mp.polylog(n, z))

Li5v = {'b': Lin(5,b), '-b': Lin(5,-b), '1/2': Lin(5,mp.mpf('0.5')), '1/sqrt2': Lin(5,1/r2),
        '2-sqrt2': Lin(5,2-r2), 'p': Li_real(5, pval)}
Li4v = {'b': Lin(4,b), '-b': Lin(4,-b), '1/2': Lin(4,mp.mpf('0.5')), '1/sqrt2': Lin(4,1/r2),
        '2-sqrt2': Lin(4,2-r2), 'p': Li_real(4, pval)}
Li3v = {'b': Lin(3,b), '1/2': Lin(3,mp.mpf('0.5')), '1/sqrt2': Lin(3,1/r2)}
Li2v = {'b': Lin(2,b)}
# cross-check real branch Li_2 at p>1 vs real integral int_0^1 log|1-p t|/t dt  (= -Li_2^real(p))
mp.mp.dps = 60
chk = -mp.quad(lambda t: mp.log(abs(1-pval*t))/t, [0, 1/pval, 1], method='tanh-sinh')
print("Re(polylog(2,p)) vs real integral:", mp.nstr(Li_real(2,pval)-chk, 4))
# and Li_4 real branch via nested real quad: G(0,0,0,a;1) = int int int log|1-r a|/(t s) dr ds dt, a=1/p
av = 1/pval
nested = mp.quad(lambda t: mp.quad(lambda s: mp.quad(lambda r: mp.log(abs(1-r*av)), [0,s])/(s), [0,t])/(t), [0,1], method='tanh-sinh')
print("Re(polylog(4,p)) vs nested real quad:", mp.nstr(Li_real(4,pval)+nested, 4))
mp.mp.dps = 250

# Li32 atoms: double series sum_{n>m>=1} x^n y^m/(n^3 m^2)
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
