import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, log, nstr, matrix, lu_solve
mp.dps = 80
data = json.load(open('ckpt_lame_p1.json'))
def extract(taus, Ps, P0, lead, nterms):
    n = nterms
    Q = [(P-P0)/t**lead for P,t in zip(Ps,taus)]
    A = matrix(n,n); b = matrix(n,1)
    for i,t in enumerate(taus):
        for j in range(n): A[i,j] = t**j
        b[i] = Q[i]
    return lu_solve(A,b)
for lam in ['2','3','0.5']:
    epss = [mpf(s) for s in data[f'lam{lam}_pinf_epss']]
    Ps = [mpf(s) for s in data[f'lam{lam}_pinf_Ps']]
    P0 = mpf(8)
    c = extract(epss[2:9], Ps[2:9], P0, 1, 7)
    print(f"--- lam={lam} p->inf (eps=1/p):")
    for j,ci in enumerate(c):
        print(f"  d{j+1} = {nstr(ci,30)}")
