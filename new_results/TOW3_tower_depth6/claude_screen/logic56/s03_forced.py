# s03: FIT-FREE forced recursion, non-adjacent, symbolic (d,m):
#   F(delta,0) := a_delta (rho/e00, closed forms; delta<=4 TOW2 cache, 5,6 Bell)
#   F(delta,s+1) := Lambda(delta,s) - F(delta,s)      (s = 0..floor(delta/2))
# Closing conditions: F(delta, floor(delta/2)+1) == 0.  Compare forced F with the lib's fitted/zigzag forms.
import sys, time
sys.dont_write_bytecode = True
import sympy as sp
import tow2lib as T2, tow3lib as T3
d, m = T3.d, T3.m
I = sp.I
tau1 = T2.tau1f; tau2 = T2.tau2f
t0 = time.time()
Fz = {}
def F(dp, s):
    """forced family; zero above the closing level (those zeros are checked below)."""
    if s > dp//2: return sp.Integer(0)
    return Fz[(dp, s)]
def Lam_own(delta, s):
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        fA = F(dp, s).subs({d: d+1, m: m+1}, simultaneous=True)
        fB = F(dp, s).subs({d: d+1, m: m}, simultaneous=True)
        tot += I*(2*d-1)/(2*d-dp-1) * (tau1(delta-dp)*fA - tau2(delta-dp)*fB)
    return tot
closing = {}
for delta in range(0, 7):
    Fz[(delta, 0)] = sp.cancel(T3.arho(delta))
    for s in range(0, delta//2+1):
        nxt = sp.cancel(sp.together(Lam_own(delta, s) - F(delta, s)))
        if s+1 <= delta//2:
            Fz[(delta, s+1)] = nxt
        else:
            closing[delta] = (nxt == 0)
            print(f'  closing condition delta={delta}: F(delta,{s+1}) == 0 : {nxt==0}', flush=True)
    print(f'delta={delta} done ({round(time.time()-t0,1)}s)', flush=True)
print('--- compare forced vs library (fitted / zigzag / TOW2) forms ---')
allok = True
for delta in range(0, 7):
    for s in range(0, delta//2+1):
        diff = sp.cancel(sp.together(Fz[(delta, s)] - T3.e(delta, s)))
        ok = (diff == 0)
        allok &= ok
        tag = 'fitted' if (delta, s) in [(4,1),(5,1),(6,1),(6,2)] else ('top' if s == delta//2 and s>0 else ('rho' if s==0 else 'other'))
        print(f'  forced F({delta},{s}) == lib e({delta},{s})  [{tag}] : {ok}')
print('ALL closing conditions:', all(closing.values()), '  ALL forced==lib:', allok, round(time.time()-t0,1), 's')
import pickle
pickle.dump({k: str(v) for k, v in Fz.items()}, open('forced_F.pkl', 'wb'))
