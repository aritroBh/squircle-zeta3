import mpmath
from mpmath import mp, mpf, sqrt, pi, zeta, gamma, log, cos, sin, exp
mp.dps = 75

def P_squircle_tau(tau):
    """Perimeter of |x|^p+|y|^p=1 with tau = 1-1/p (diamond end tau->0+)."""
    p = 1/(1-tau)
    q = p/(p-1)   # = 1/tau
    f = lambda th: (cos(th)**q + sin(th)**q)**(1/q)
    return 8*mpmath.quad(f, [0, pi/4])

# sample at geometrically spaced tau
taus = [mpf('0.06')*mpf('0.72')**j for j in range(18)]
Ps = [P_squircle_tau(t) for t in taus]
for t,P in zip(taus[:6],Ps[:6]):
    print(t, P)
import json
json.dump([str(t) for t in taus], open('ckpt_squircle_taus.json','w'))
json.dump([str(P) for P in Ps], open('ckpt_squircle_Ps.json','w'))
print("saved")
