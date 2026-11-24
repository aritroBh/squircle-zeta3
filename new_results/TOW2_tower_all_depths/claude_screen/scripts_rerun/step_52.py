import sys, os, time
work = '<project>/research/08_theory/TOW2_tower_all_depths/claude_screen/scripts_rerun'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
d, m = T.d, T.m
i = T.i

def e_closed(delta, s, dd, mm):
    """e(delta,s)/e00 at (d=dd, m=mm), off-diagonal closed forms."""
    dd_, mm_ = T.d, T.m
    sub = {dd_: dd, mm_: mm}
    if s == 0: return T.arho(delta).subs(sub)
    if delta == 2 and s == 1: return (-1/(2*dd-3))
    if delta == 3 and s == 1: return i*(dd-3)/(2*dd-3)
    if delta == 4 and s == 1:
        g = dd**3 - 9*dd**2 + sp.Rational(113,4)*dd - i*dd/2 + i*mm - 29
        return g/((2*dd-3)*(2*dd-5))
    if delta == 4 and s == 2: return sp.Rational(2)/((2*dd-3)*(2*dd-5))
    raise ValueError((delta,s))

# (4,1) slope identity, non-adjacent:
lhs = e_closed(4,1,d,m) + e_closed(4,2,d,m)
rhs = 0
for dp in range(0,5):
    if dp < 2: continue   # parents' e(dp,1)=0 for dp<2
    cA = T.tau1(4-dp)*(m+1)/(2*d-dp-1); cB = T.tau2(4-dp)*(d-m+1)/(2*d-dp-1)
    rhs += cA*e_closed(dp,1,d+1,m+1)*T.rA0() + cB*e_closed(dp,1,d+1,m)*T.rB0()
t0 = time.time()
rem = sp.cancel(sp.together(lhs - rhs))
print('E4-slope (4,1) non-adjacent remainder:', rem, ' (', round(time.time()-t0,1),'s )')
