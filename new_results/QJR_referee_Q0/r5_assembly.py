# Referee QJR: exact assembly 6(P1-P2a-P2c) - Q0 (theorem) = 0, using r4 (P1 reduced, own) and own P2a/P2c reductions.
from rcommon import *
import json
P1=sp.sympify(json.load(open('r4_P1_reduced.json')))
P2a=r2*(3*Lam-M/4-lam*L/2+L*l3+L**2/2-5*pi**2/24)      # reproduced exactly in r2_P2.log
P2c=r2*(Lam-M/8+L*l3/2+pi**2/96-pi*alpha/4)              # reproduced exactly in r2_P2.log
Q0=(-pi**2/8+r2/16*pi**2+3*r2/2*pi*alpha-3*lam**2-9*r2*L*l3+(6+6*r2)*L**2-6*r2*Lam+9*r2/4*M)
d=sp.expand(6*(P1-P2a-P2c)-Q0)
print('6(P1-P2a-P2c) - Q0 =',d)
mp.mp.dps=110
v=numval(Q0.subs({Lam:Li2(b),M:Li2(-b**4)}))
print('Q0 =',mp.nstr(v,100))
