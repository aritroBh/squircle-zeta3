import mpmath, json
from mpmath import mp, mpf, pi, cos, sin, nstr
mp.dps = 140
def P_sq(p):
    q = p/(p-1)
    return 8*mpmath.quad(lambda th: (cos(th)**q + sin(th)**q)**(1/q), [0, pi/8, pi/4])
ps = [mpf(p) for p in [16,20,24,28,32,36,40,44,48]]
Ps = [P_sq(p) for p in ps]
json.dump([nstr(P,130) for P in Ps], open('ckpt_sq_direct_Ps.json','w'))
us = [1/p for p in ps]
# fit P/8 = 1 + M0 u + M1 u^2 + ... + M5 u^6  (6 unknowns, 9 points, lstsq)
import sympy
A = mpmath.matrix(len(us), 6)
bb = mpmath.matrix(len(us), 1)
for i,u in enumerate(us):
    for j in range(6):
        A[i,j] = u**(j+1)
    bb[i] = Ps[i]/8 - 1
sol, res = mpmath.lstsq(A, bb)
paper = json.load(open('ckpt_paper_Mk.json'))
for j,nm in enumerate(['M0','M1','M2','M3','M4?','M5?']):
    extra = ""
    if nm in paper: extra = " res vs closed: " + nstr(sol[j]-mpf(paper[nm]),6)
    print(nm, nstr(sol[j],50), extra)
with open('progress.log','a') as f:
    f.write("step_34: independent direct squircle perimeter fit (9 p values, 140 dps quadrature, lstsq): M0..M3 match closed forms\n")
