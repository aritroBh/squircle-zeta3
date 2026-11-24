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
out = {}
for lam in ['2','3','0.5']:
    epss = [mpf(s) for s in data[f'lam{lam}_pinf2_epss']]
    Ps = [mpf(s) for s in data[f'lam{lam}_pinf2_Ps']]
    c = extract(epss, Ps, mpf(8), 1, 6)
    out[lam] = [nstr(ci,60) for ci in c]
    print(f"--- lam={lam} p->inf:")
    for j,ci in enumerate(c):
        print(f"  d{j+1} = {nstr(ci,35)}")
json.dump(out, open('ckpt_lame_pinf_coeffs.json','w'))
