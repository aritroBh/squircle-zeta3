helpers = '''

def polyval(key, delta, k, diag_shift=None):
    """tower-law polynomial for key at depth delta, evaluated at level k.
       Fit from levels k0..k0+deg where k0 = d (off-diag) or d+1 (diagonal)."""
    d = key[0]+key[1]
    if diag_shift is None:
        diag_shift = (key[0]==key[1])
    k0 = d + (1 if diag_shift else 0)
    degf = delta//2
    xs = list(range(k0, k0+degf+1))
    ys = [u(key[0],key[1],x,delta) for x in xs]
    out = ZERO
    for a,xa in enumerate(xs):
        La = Fr(1)
        for b,xb in enumerate(xs):
            if b!=a: La *= Fr(k-xb, xa-xb)
        out = out + ys[a].scale(La)
    return out

def feed_diag_into_child(p, k, delta, Udiag):
    i = p; j = p-1
    return Lambda_ij(i,j,k,delta, lambda kk,dp: Udiag(dp) if kk==(p,p) else ZERO)
'''
open('towlib.py','a').write(helpers)
import importlib, towlib
importlib.reload(towlib)
from towlib import *
for p in (2,3,4):
    maxd = min(4*p-4, 2*(2*p-1)-4, 2*(14-2*p-1))
    print(f'--- p={p} (max depth tested {maxd}) ---')
    Deltas = {}
    for dp in range(0, maxd+1):
        D = rho((p,p),dp) - polyval((p,p),dp,2*p)
        Deltas[dp] = D
        kind = 'real' if D.is_real() else ('pure-i' if (D.p0==0 and D.p1==0) else 'COMPLEX')
        print(f'  Delta({dp}) = {D}  [{kind}]')
    for delta in range(0, maxd+1):
        feed = feed_diag_into_child(p, 2*p, delta, lambda dp: Deltas.get(dp, ZERO))
        print(f'  feed into child at delta={delta}: zero={feed.is_zero()}  {"" if feed.is_zero() else feed}')
