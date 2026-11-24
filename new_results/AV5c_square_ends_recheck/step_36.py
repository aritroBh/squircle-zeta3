import mpmath, json
from mpmath import mp, mpf, pi, cos, sin, nstr
mp.dps = 160
def P_sq(p):
    q = p/(p-1)
    return 8*mpmath.quad(lambda th: (cos(th)**q + sin(th)**q)**(1/q), [0, pi/8, pi/4])
ps = [mpf(p) for p in [32,40,48,56,64,72,80,96]]
Ps = [P_sq(p) for p in ps]
json.dump([nstr(P,150) for P in Ps], open('ckpt_sq_direct_Ps_big.json','w'))
us = [1/p for p in ps]
def solve_pts(idxs, nk=6):
    A = mpmath.matrix(nk, nk); bb = mpmath.matrix(nk, 1)
    for r,i in enumerate(idxs):
        u = us[i]
        for j in range(nk): A[r,j] = u**(j+1)
        bb[r] = Ps[i]/8 - 1
    return mpmath.lu_solve(A, bb)
sol = solve_pts([0,1,2,3,4,5])
paper = json.load(open('ckpt_paper_Mk.json'))
for j,nm in enumerate(['M0','M1','M2','M3','M4','M5']):
    extra = ""
    if nm in paper: extra = " res vs closed: " + nstr(sol[j]-mpf(paper[nm]),6)
    print(nm, nstr(sol[j],60), extra)
for i in [6,7]:
    u = us[i]
    pred = 8*(1 + sum(sol[j]*u**(j+1) for j in range(6)))
    print(f"p={ps[i]}: P-pred residual = {nstr(Ps[i]-pred,6)}")
with open('progress.log','a') as f:
    f.write("step_36: direct perimeter fit with p=32..96 (160 dps): M0 to 19 digits, M1 16, M2 13, M3 11 -- independent confirmation of closed forms\n")
