import mpmath, pickle, json
from mpmath import mp, mpf, sqrt, pi, cos, sin, matrix, lu_solve, nstr
mp.dps = 90
def P_sq_support(tau):
    q = 1/tau
    f = lambda th: (cos(th)**q + sin(th)**q)**(1/q)
    return 8*mpmath.quad(f, [0, pi/8, pi/4])
taus = [mpf('0.02')/mpf('1.85')**j for j in range(12)]
Ps = [P_sq_support(t) for t in taus]
json.dump([nstr(t,30) for t in taus], open('ckpt_sq2_taus.json','w'))
json.dump([nstr(P,80) for P in Ps], open('ckpt_sq2_Ps.json','w'))
print("computed", [nstr(t,8) for t in taus])
print(nstr(Ps[-1],60))
