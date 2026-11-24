import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, cos, sin, matrix, lu_solve
mp.dps = 75

def P_sq_support(tau):
    p = 1/(1-tau); q = 1/tau
    f = lambda th: (cos(th)**q + sin(th)**q)**(1/q)
    return 8*mpmath.quad(f, [0, pi/4])

def extract(taus, Ps, P0, lead, nterms, verbose=True):
    """fit (P-P0)/tau^lead = sum_{j=0}^{nterms-1} d_j tau^j by interpolation"""
    n = nterms
    Q = [(P-P0)/t**lead for P,t in zip(Ps,taus)]
    A = matrix(n,n); b = matrix(n,1)
    for i,t in enumerate(taus):
        for j in range(n): A[i,j] = t**j
        b[i] = Q[i]
    c = lu_solve(A,b)
    # residual check on all data
    r = matrix(n,1)
    for i in range(n):
        pred = sum(c[j]*taus[i]**j for j in range(n))
        r[i] = pred - Q[i]
    return c, r

taus6 = [mpf('0.032')/2**j for j in range(6)]
Ps6 = [P_sq_support(t) for t in taus6]
P0 = 4*sqrt(2)
c, r = extract(taus6, Ps6, P0, 2, 6)
for j,ci in enumerate(c):
    print(f"c{j+2} = {mp.nstr(ci,45)}")
print("max resid:", mp.nstr(max(abs(x) for x in r),4))
print("pred c2 =", mp.nstr(sqrt(2)*zeta(2),45))
print("pred c3 =", mp.nstr(sqrt(2)*zeta(3),45))
json.dump([str(t) for t in taus6], open('ckpt_sq_taus6.json','w'))
json.dump([mp.nstr(P,70) for P in Ps6], open('ckpt_sq_Ps6.json','w'))
