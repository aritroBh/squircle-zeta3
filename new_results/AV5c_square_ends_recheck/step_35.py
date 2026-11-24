import mpmath, json
from mpmath import mp, mpf, pi, cos, sin, nstr
mp.dps = 140
Ps = [mpf(s) for s in json.load(open('ckpt_sq_direct_Ps.json'))]
ps = [mpf(p) for p in [16,20,24,28,32,36,40,44,48]]
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
    print(nm, nstr(sol[j],50), extra)
# validation: predict P at remaining points
for i in [6,7,8]:
    u = us[i]
    pred = 8*(1 + sum(sol[j]*u**(j+1) for j in range(6)))
    print(f"p={ps[i]}: P-pred residual = {nstr(Ps[i]-pred,6)}")
with open('progress.log','a') as f:
    f.write("step_35: independent direct squircle perimeter fit (140 dps quadrature, 6x6 solve + 3 validation points): M0..M3 match closed forms\n")
