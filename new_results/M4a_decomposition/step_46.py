import sys, pickle
sys.dont_write_bytecode = True
sys.path.insert(0,'.')
import mpmath as mp
from flint import acb, ctx
import dilog_exact as DE

mp.mp.dps = 1250
ctx.prec = 4200
tails = pickle.load(open('m3_tails_exact.pkl','rb'))
W3p = tails['W3p']
print('structure:', W3p)

r2 = mp.sqrt(2)
# mpmath evaluation
v_mp = W3p.num()
# flint evaluation of the same expression
base_fl = {'lam': acb(2).log(), 'L': (acb(2).sqrt()+1).log(), 'pi2': acb.pi()**2, 'z3': acb(3).zeta()}
tot = acb(0)
for k, c in W3p.d.items():
    x = acb(str(c.p))/1 + acb(2).sqrt()*acb(str(c.q))
    for a in k:
        if isinstance(a, str):
            x *= base_fl[a]
        else:
            z = DE.float_exact(a[1])
            x *= acb(mp.nstr(z, 1200)).polylog(acb(2 if a[0]=='Li2' else 3))
    tot += x
v_fl = mp.mpf(tot.real.mid().str(1200, radius=False, more=True))
print('mp  =', mp.nstr(v_mp, 100))
print('fl  =', mp.nstr(v_fl, 100))
print('mp - fl =', mp.nstr(v_mp - v_fl, 5))
print('flint rad:', tot.real.rad().str(5))
