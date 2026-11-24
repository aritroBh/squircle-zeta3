import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, log, nstr, matrix, lu_solve, asinh
mp.dps = 80
data = json.load(open('ckpt_lame_p1.json'))
def P1_lam(lam):  # perimeter at p=1 exactly
    if lam=='2' or lam=='0.5': return 2*(sqrt(5)+asinh(2)/2)
    # lam=3: numeric
    a=2; b=mpf(2)/3
    def f(t):
        dx = -a*cos(t)**(a-1)*sin(t)
        dy = b*sin(t)**(b-1)*cos(t)
        return sqrt(dx**2+dy**2)
    return 4*mpmath.quad(f, [0, pi/4, pi/2])
def extract(taus, Ps, P0, lead, nterms):
    n = nterms
    Q = [(P-P0)/t**lead for P,t in zip(Ps,taus)]
    A = matrix(n,n); b = matrix(n,1)
    for i,t in enumerate(taus):
        for j in range(n): A[i,j] = t**j
        b[i] = Q[i]
    return lu_solve(A,b)
for lam in ['2','3','0.5']:
    taus = [mpf(s) for s in data[f'lam{lam}_p1_taus']]
    Ps = [mpf(s) for s in data[f'lam{lam}_p1_Ps']]
    P0 = P1_lam(lam)
    # try lead=1 (analytic) with 6 terms
    c = extract(taus[2:9], Ps[2:9], P0, 1, 7)
    print(f"--- lam={lam}, P0={nstr(P0,30)}")
    for j,ci in enumerate(c):
        print(f"  c{j+1} = {nstr(ci,30)}")
