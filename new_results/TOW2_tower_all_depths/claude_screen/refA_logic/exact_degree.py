import sys; sys.path.insert(0,'.')
import sympy as sp
import tow2lib as T
i = sp.I; p = sp.Symbol('p')
def Fq(delta, s, dd, mm):
    if s > delta//2: return 0
    if s == 0: return T.arho_simple(delta).subs({T.d: dd, T.m: mm}, simultaneous=True) if delta else 1
    if (delta,s)==(2,1): return -1/(2*dd-3)
    if (delta,s)==(3,1): return i*(dd-3)/(2*dd-3)
    if (delta,s)==(4,1): return (dd**3-9*dd**2+sp.Rational(113,4)*dd-i*dd/2+i*mm-29)/((2*dd-3)*(2*dd-5))
    if (delta,s)==(4,2): return sp.Rational(2)/((2*dd-3)*(2*dd-5))
def ehat(delta, s):
    return sp.cancel(sp.together((p+1)*sum(T.tau1f(delta-dp)*Fq(dp, s, 2*p+1, p+1)*(4*p-1)*(4*p-3)/((p+1)*p)/(4*p-dp-1) for dp in range(delta+1))))
for (delta, s) in [(2,1),(3,1),(4,2)]:
    e = sp.factor(ehat(delta, s))
    print(f'diag top e-hat({delta},{s}) (units e00_(p,p-1)) =', e, '| roots:', sp.solve(sp.numer(sp.together(e)), p))
for (delta, s) in [(2,1),(3,1),(4,2)]:
    print('NA top F(',delta,s,') =', sp.factor(Fq(delta, s, sp.Symbol('d'), 0)), '(nonzero for d>=4)')
