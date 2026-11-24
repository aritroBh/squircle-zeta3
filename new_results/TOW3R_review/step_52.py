work = '<project>/research/08_theory/TOW3R_review/'
src = '''import sys, os, time, pickle
from fractions import Fraction as Fr
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work+'/tow3_copy')
import sympy as sp
import tow3lib as T3
d, m = T3.d, T3.m; i = sp.I

# (a) e_fits.pkl vs tow3lib closed forms
ef = pickle.load(open('e_fits.pkl','rb'))
def fit_expr(rec):
    b = rec['b']; deg = rec['deg']; mons = rec['mons']; sol = rec['sol']
    num = sum((re_+im_*i)*d**di*m**mi for (di,mi),(re_,im_) in zip(mons,sol))
    return sp.cancel(num / T3.Pi(b))
for key in [(4,1),(5,1),(6,2)]:
    fe = fit_expr(ef[key])
    lib = T3._EF[key]
    print(f'e{key}: fit pickle == tow3lib closed form: {sp.expand(fe-lib)==0}', flush=True)
e61 = sp.sympify(pickle.load(open('e61_closed.pkl','rb')), locals={'d':d,'m':m,'I':i})
print('e(6,1): e61_closed == tow3lib:', sp.expand(sp.cancel(e61)-T3.e61)==0, flush=True)
# also writeup e(5,1), e(6,2) denominators: 12 Pi_2, 2 Pi_3 -- compare with fits' Pi_b
print('fit denominators: e(4,1): Pi_%d, e(5,1): Pi_%d, e(6,2): Pi_%d' % (ef[(4,1)]['b'], ef[(5,1)]['b'], ef[(6,2)]['b']))
# NOTE: writeup says e(5,1) = N/(12 Pi_2); fit has b=2 => denominator Pi_2; check ratio
f51 = fit_expr(ef[(5,1)])
print('e(5,1) fit numerator/denominator vs lib (lib = N51/(12 Pi_2)): equal already checked above')

# (b) extract ACTUAL Newton coefficients from tower data and compare, all keys, depths 5,6
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def gdiv(x,y):
    den = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/den, (x[1]*y[0]-x[0]*y[1])/den)
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gadd(x,y): return (x[0]+y[0],x[1]+y[1])
def gsub(x,y): return (x[0]-y[0],x[1]-y[1])
def qsi_div(x, y):
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    den = (gmul(C,C)[0]-2*gmul(D,D)[0], gmul(C,C)[1]-2*gmul(D,D)[1])
    numA = (gmul(A,C)[0]-2*gmul(B,D)[0], gmul(A,C)[1]-2*gmul(B,D)[1])
    numB = (gmul(B,C)[0]-gmul(A,D)[0], gmul(B,C)[1]-gmul(A,D)[1])
    return (gdiv(numA, den), gdiv(numB, den))
def qsi_gmul(x, c):  # x QSI times Gaussian rational c -> QSI
    A=(x[0],x[2]); B=(x[1],x[3])
    # (A + B s2) c
    return (gmul(A,c)[0], gmul(B,c)[0], gmul(A,c)[1], gmul(B,c)[1])
from math import comb
t0=time.time()
res = {}
for delta in [5,6]:
    smax = delta//2
    bad = []; tested = 0
    for key in H:
        if key=='_meta': continue
        p,q=key; dd=p+q
        if q==p:  # diagonal: Newton basis is shifted; skip here (checked separately)
            continue
        if delta not in H.get(key,{}) or 2*dd-4 < delta: continue
        e00 = H[key][0][dd]
        # actual e(delta,s): u(delta,k) = sum_s binom(k-dd,s) e00*e(s); invert using k=dd..dd+smax
        vals = []
        okk = True
        for kk in range(dd, dd+smax+1):
            if kk not in H[key][delta]: okk=False; break
            vals.append(qsi_div(H[key][delta][kk], e00))
        if not okk: continue
        # invert: v_t = sum_{s<=t} binom(t,s) E_s  => E_t = v_t - sum_{s<t} binom(t,s) E_s
        E = []
        for tt in range(smax+1):
            acc = vals[tt]
            for s in range(tt):
                c = (Fr(comb(tt,s)), Fr(0))
                term = qsi_gmul(E[s], c)
                acc = (gsub(acc[0],term[:2:1]+(term[1],)) if False else acc)
                acc = ((acc[0][0]-term[0][0], acc[0][1]-term[0][1]), (acc[1][0]-term[1][0], acc[1][1]-term[1][1]))
            E.append(acc)
        for s in range(smax+1):
            want = sp.expand(T3.e(delta,s).subs({d:dd, m:p}))
            wr, wi = sp.Rational(sp.re(want)), sp.Rational(sp.im(want))
            got = E[s]
            tested += 1
            if got[1]!=(Fr(0),Fr(0)) or got[0]!=(wr,wi):
                bad.append((key,s,got,(wr,wi)))
    res[delta]=(tested,len(bad))
    print(f'depth {delta}: off-diagonal ACTUAL Newton coefficients vs closed forms: tested {tested}, bad {len(bad)} {bad[:3]}  ({time.time()-t0:.0f}s)', flush=True)
'''
open(work+'step_12.py','w').write(src)
import subprocess, datetime
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,work+'step_12.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout, r.stderr[-2000:])
open(work+'step_12.out','w').write(r.stdout+r.stderr)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_12 fits==lib; actual Newton coefs from data vs closed forms depths 5,6\n")
