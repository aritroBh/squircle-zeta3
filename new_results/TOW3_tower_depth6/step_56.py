import sys, os, time, pickle
from fractions import Fraction as Fr
from math import comb
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
os.nice(10)
import sympy as sp
import tow3lib as T3
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def u(i,j,k,delta):
    return tuple(H[(i,j)][delta][k]) if (i,j) in H and delta in H[(i,j)] and k in H[(i,j)][delta] else (Fr(0),)*4
def qsi_scale_g(x, g):  # multiply QSI x by Gaussian rational g=(re,im)
    A=(x[0],x[2]); B=(x[1],x[3])
    def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
    A2 = gmul(A,g); B2 = gmul(B,g)
    return (A2[0], B2[0], A2[1], B2[1])
def qsi_add(x,y): return tuple(a+b for a,b in zip(x,y))
def g_of_expr(expr, dd, mm):
    v = sp.expand(expr.subs({T3.d: dd, T3.m: mm}))
    return (sp.Rational(sp.re(v)), sp.Rational(sp.im(v)))

t0 = time.time()
bad = []; tested = 0
for key in H:
    if key == '_meta': continue
    p, q = key; dd = p+q
    if 2*dd-4 < 0: continue
    maxdelta = min(6, 2*dd-4)
    e00 = u(p,q,dd,0)
    # precompute e-ratio Gaussian values
    erat = {}
    for delta in range(0, maxdelta+1):
        for s in range(0, delta//2+1):
            erat[(delta,s)] = g_of_expr(T3.e(delta,s), dd, p)
    if p == q:
        # ehat ratios
        ehat = {}
        for delta in range(0, maxdelta+1):
            for s in range(0, delta//2+1):
                tot = (Fr(0),Fr(0))
                for dpp in range(0, delta+1):
                    g = g_of_expr(T3.e(dpp,s), 2*dd+1, dd+1)
                    t1 = T3.tau1(delta-dpp)
                    t1g = (sp.Rational(sp.re(t1)), sp.Rational(sp.im(t1)))
                    rB0d = (Fr(0), Fr(2*(4*dd-1), dd+1))
                    c = (sp.Rational(sp.re(t1)), sp.Rational(sp.im(t1)))
                    # g * tau1 * rB0d / (4p-dpp-1) * (p+1)
                    def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
                    term = gmul(gmul(g, t1g), rB0d)
                    term = (term[0]*Fr(dd+1,4*dd-dpp-1), term[1]*Fr(dd+1,4*dd-dpp-1))
                    tot = (tot[0]+term[0], tot[1]+term[1])
                ehat[(delta,s)] = tot
    for k in range(dd, min(dd+8, 52)+1):
        if p == q and k < 2*p+1: continue
        for delta in range(0, maxdelta+1):
            U = (Fr(0),)*4
            if p == q:
                for s in range(0, delta//2+1):
                    U = qsi_add(U, qsi_scale_g(e00, (ehat[(delta,s)][0]*comb(k-2*p-1,s), ehat[(delta,s)][1]*comb(k-2*p-1,s))))
            else:
                for s in range(0, delta//2+1):
                    U = qsi_add(U, qsi_scale_g(e00, (erat[(delta,s)][0]*comb(k-dd,s), erat[(delta,s)][1]*comb(k-dd,s))))
            actual = u(p,q,k,delta)
            tested += 1
            if U != actual:
                bad.append((key,k,delta))
print(f'tower law delta<=6 vs data: {tested} cells, {len(bad)} bad  ({round(time.time()-t0,1)}s)')
if bad[:5]: print(bad[:5])
with open(f'{work}/progress.log','a') as f:
    f.write(f'step_56: full tower law delta<=6 vs TOWN1 data: {tested} cells, {len(bad)} bad\n')
