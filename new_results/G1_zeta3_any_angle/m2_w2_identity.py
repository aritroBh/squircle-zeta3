"""m2_w2_identity.py -- the weight-2 part.  From m2_structure.py, W2 = -int_1^{t1} log g dt/(t+1).
Explicit antiderivative (m2_split_check.W2_closed) and the reduction
   W2 = -I_1 - I_{-1} + I_ca + log2 log(1+s)  =  (1/2) Li2(s^2) + 2 log s log c  =  zeta(2)/2 - Li2(c^2)/2
(Landen: Li2(-s/(1-s)) = -Li2(s) - log^2(1-s)/2, Li2(-s^2/c^2) = -Li2(s^2) - log^2(c^2)/2,
 duplication Li2(s)+Li2(-s) = Li2(s^2)/2, Euler: Li2(s^2)+Li2(c^2) = z2 - log s^2 log c^2).
The log bookkeeping is done symbolically below (formal logs A=log2, P=log(1+s), M=log(1-s), S=log s,
log c^2 = P+M), then everything is checked numerically at the 9 angles."""
import sympy as sp
import mpmath as mp
A, P, M, S = sp.symbols('A P M S')
logc2 = P + M
I1_logs = (A + S)*P                                  # log(2s) log(1+s)
Im1 = ((A + P)**2 - A**2)/2
Ica_logs = (A + S + P)*(P - logc2) - M**2/2 - (A + 2*S)*(-logc2) + logc2**2/2
logs = sp.expand(-I1_logs - Im1 + Ica_logs + A*P)
print('log part of W2 =', sp.factor(logs), '  (expected S*(P+M) = log s * log c^2)')
from m2_split_check import W2_closed, ANG
mp.mp.dps = 120
for tag, a in ANG.items():
    c = mp.cos(a/2)
    print(f'{tag:>7}: |W2_closed - (z2/2 - Li2(c^2)/2)| = {mp.nstr(abs(W2_closed(a) - (mp.zeta(2) - mp.polylog(2, c*c))/2), 3)}')
