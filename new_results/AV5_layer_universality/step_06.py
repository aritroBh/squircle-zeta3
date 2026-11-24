import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, log, cos, sin, exp, matrix, lu_solve
mp.dps = 75
taus = [mpf(s) for s in json.load(open('ckpt_squircle_taus.json'))]
Ps   = [mpf(s) for s in json.load(open('ckpt_squircle_Ps.json'))]

def fit_coeffs(taus, vals, P0, powers):
    """Least-squares fit vals - P0 = sum c_k tau^k over given powers."""
    n = len(powers)
    # normal equations with column scaling
    scales = [max(abs(t)**k for t in taus) for k in powers]
    A = matrix(len(taus), n)
    for i,t in enumerate(taus):
        for j,k in enumerate(powers):
            A[i,j] = t**k / scales[j]
    b = matrix(len(taus),1)
    for i,v in enumerate(vals): b[i,0] = v - P0
    # solve A^T A c = A^T b
    AtA = matrix(n,n); Atb = matrix(n,1)
    for i in range(n):
        for j in range(n):
            AtA[i,j] = sum(A[r,i]*A[r,j] for r in range(len(taus)))
        Atb[i,0] = sum(A[r,i]*b[r,0] for r in range(len(taus)))
    c = lu_solve(AtA, Atb)
    return [c[j,0]*scales[j] for j in range(n)]

P0 = 4*sqrt(2)
powers = [2,3,4,5,6,7]
cs = fit_coeffs(taus, Ps, P0, powers)
for k,c in zip(powers,cs):
    print(f"c{k} = {mp.nstr(c,50)}")
print()
print("predictions: c2 = 2sqrt2 zeta2 =", mp.nstr(2*sqrt(2)*zeta(2),50))
print("             c3 = 2sqrt2 zeta3 =", mp.nstr(2*sqrt(2)*zeta(3),50))
print("ratios: c2/(2sqrt2)=", mp.nstr(cs[0]/(2*sqrt(2)),40))
print("        c3/(2sqrt2)=", mp.nstr(cs[1]/(2*sqrt(2)),40))
json.dump([mp.nstr(c,60) for c in cs], open('ckpt_squircle_coeffs.json','w'))
