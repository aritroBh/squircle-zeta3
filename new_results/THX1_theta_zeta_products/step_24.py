
import sympy as sp
import mpmath as mp
mp.mp.dps = 60
import pickle
from fractions import Fraction as Fr

# Profile: P = sqrt((1+y)/2) * (1 + t*t^w)^{1/(2(1+w))} - 1,  t=(1-y)/(1+y), w = tau/(1-tau)
# log P1 = log(1 + t e^{w L}) / (2(1+w)), L = log t
w, L, M = sp.symbols('w L M')   # L = log t, M = log(1+t) = log(2/(1+y))
# series of (1/(2(1+w))) * log(1 + t e^{wL}) in w; write log(1+t e^{wL}) = log(1+t) + log(1 + t(e^{wL}-1)/(1+t))
r = sp.symbols('r')  # r = t/(1+t) = (1-y)/2
t = sp.symbols('t')
N = 8
ew = sum((w*L)**k/sp.factorial(k) for k in range(N+1))
inner = sp.series(sp.log(1 + r*(ew-1)), w, 0, N).removeO().expand()
logP1 = sp.expand(inner + M) * sp.series(1/(2*(1+w)), w, 0, N).removeO().expand()
logP1 = sp.expand(logP1)
def truncW(expr, n):
    expr = sp.expand(expr)
    return sum(expr.coeff(w, i)*w**i for i in range(n+1))
logP1 = truncW(logP1, N)
# P1 = exp(logP1) - 1
P1 = sp.exp(logP1) - 1
P1 = truncW(sp.expand(P1, expt=False) if False else sp.series(sp.exp(logP1), w, 0, N).removeO() - 1, N)
# substitute r = t/(1+t) is already the symbol r; note we treat M, L, r independent: M=log(1+t), L=log t, r=t/(1+t)
print("P1 (in w): order 1:", sp.expand(P1).coeff(w,1))
print("order 2:", sp.expand(P1).coeff(w,2))
# w = tau/(1-tau): re-expand
tau = sp.symbols('tau')
Ptau = sp.series(P1.subs(w, tau/(1-tau)), tau, 0, 7).removeO().expand()
for n in range(1,6):
    print("tau^%d coeff:"%n, sp.collect(sp.expand(Ptau).coeff(tau,n), [L,M]))
