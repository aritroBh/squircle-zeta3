import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, log, cos, sin, nstr, matrix, lu_solve, catalan
mp.dps = 100
def P_sq_eps(eps):
    q = 1+eps
    f = lambda th: (cos(th)**q + sin(th)**q)**(1/q)
    return 8*mpmath.quad(f, [0, pi/4])
epss = [mpf('0.02')/mpf('1.8')**j for j in range(12)]
Ps = [P_sq_eps(e) for e in epss]
json.dump([nstr(e,30) for e in epss], open('ckpt_sqsq_epss.json','w'))
json.dump([nstr(P,90) for P in Ps], open('ckpt_sqsq_Ps.json','w'))
print(Ps[0], Ps[-1])
# extract d1..d8 from (P-8)/eps
def extract(taus, Ps, P0, nterms):
    n = nterms
    Q = [(P-P0)/t for P,t in zip(Ps,taus)]
    A = matrix(n,n); b = matrix(n,1)
    for i,t in enumerate(taus):
        for j in range(n): A[i,j] = t**j
        b[i] = Q[i]
    return lu_solve(A,b)
sel = list(range(3,11))
c = extract([epss[i] for i in sel],[Ps[i] for i in sel],8,8)
for j,ci in enumerate(c):
    print(f"d{j+1} = {nstr(ci,45)}")
