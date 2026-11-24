"""Step 3 numerics: validate the three new atom conversions and the whole Delta vector with the
project's arb level-8 evaluator (eval8, ~150 digits).  Numerical evidence only."""
import sys, pickle, time
sys.dont_write_bytecode=True
TH='<project>/research/08_theory'
sys.path.insert(0, TH+'/A3f6_level16')
from flint import ctx, acb, arb
import eval8 as ev
import mpmath as mp
ctx.prec=560; ctx.cap=1100; ev.N=1000
ev.roots=[acb(0,arb.pi()*j/4).exp() for j in range(8)]
ev.roots[0],ev.roots[2],ev.roots[4],ev.roots[6]=acb(1),acb(0,1),acb(-1),acb(0,-1)
mp.mp.dps=170
d=pickle.load(open('delta.pkl','rb'))
t0=time.monotonic()
def val(vec): 
    s=acb(0)
    for w,c in vec.items(): s+=acb(c.numerator)/c.denominator*ev.h(w)
    return s
def mpc(z): return mp.mpc(z.real.mid().str(165,radius=False), z.imag.mid().str(165,radius=False))
r2=mp.sqrt(2)
for (p,q),vec in d['new_atoms'].items():
    x=mp.mpf(p.numerator)/p.denominator+mp.mpf(q.numerator)/q.denominator*r2
    z=mpc(val(vec)); ref=mp.polylog(5,x)
    print(f'atom Li5({p}{q:+}*r2)  words={len(vec)}  |conv - mpmath| = {mp.nstr(abs(z-ref),3)}   (value {mp.nstr(ref,12)})', flush=True)
r=mpc(val(d['rat'])); s=mpc(val(d['sq']))
print('Delta rat value:', mp.nstr(r,3), ' sq value:', mp.nstr(s,3))
sc=max(abs(a) for a in d['rat'].values()); print('max |coef| rat', mp.nstr(mp.mpf(sc.numerator)/sc.denominator,4))
print(f'seconds={time.monotonic()-t0:.0f}')
