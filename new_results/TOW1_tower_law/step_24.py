from towlib import *
from fractions import Fraction as Fr

def polyval_diag(p, delta, k):
    """stable polynomial of diagonal key (p,p) at depth delta, evaluated at level k:
       fit from levels 2p+1..2p+1+floor(delta/2) using Newton/standard interpolation."""
    from fractions import Fraction as Fr
    deg = delta//2
    xs = list(range(2*p+1, 2*p+2+deg))
    ys = [u(p,p,x,delta) for x in xs]
    # Lagrange at k (componentwise on QSI 4-tuple)
    out = ZERO
    for a,xa in enumerate(xs):
        La = Fr(1)
        for b,xb in enumerate(xs):
            if b!=a: La *= Fr(k-xb, xa-xb)
        out = out + ys[a].scale(La)
    return out

def feed_diag_into_child(p, k, delta, Udiag):
    """Lambda feed into child (p,p-1) at depth delta from diagonal parent values Udiag(dp)."""
    i = p; j = p-1
    return Lambda_ij(i,j,k,delta, lambda kk,dp: Udiag(dp) if kk==(p,p) else ZERO)

for p in (2,3):
    print(f'--- p={p}: Delta(δ\') = rho_pp(δ\') - P_pp^δ\'(2p), and its feed into child ({p},{p-1}) ---')
    for dp in range(0, 2*p-3):
        Delta = rho((p,p),dp) - polyval_diag(p,dp,2*p)
        print(f'  Delta({dp}) =', Delta)
    for delta in range(0, 2*p-3):
        feed = feed_diag_into_child(p, 2*p, delta, lambda dp: rho((p,p),dp) - polyval_diag(p,dp,2*p))
        print(f'  feed at delta={delta}: {feed}  zero={feed.is_zero()}')
