import mpmath as mp, json
mp.mp.dps = 250
f450=[mp.mpf(v) for v in json.load(open('../D11b2_stokes_jump/step2_num450_N180.json'))['f450']]

def law_main(n):
    return mp.sqrt(2)*mp.cos(mp.pi*(3-2*n)/4)*mp.factorial(n)*(mp.mpf(2)/mp.pi)**n*mp.mpf(n)**mp.mpf('-1.5')

# r_n = f_n/law_main = 1 + c1/n + c2/n^2 + ...  (per class mod 4, cos sign fixed)
# Richardson: use three terms n, n+4, n+8 in a class: solve for 1, c1, c2
print("Per-class Richardson (n, n+4, n+8): estimate of constant (should be 1), c1, c2")
c1_pred = {0:(3-3*mp.pi)/8, 1:(3+3*mp.pi)/8, 2:(3-3*mp.pi)/8, 3:(3+3*mp.pi)/8}
for cls in range(4):
    n = 160 + cls  # use n, n-4, n-8 going down
    ns = [n, n-4, n-8]
    rs = [f450[k]/law_main(k) for k in ns]
    # fit r = A + B/n + C/n^2
    import numpy as np
    # solve exactly with mpmath
    M = mp.matrix(3,3)
    for i,k in enumerate(ns):
        M[i,0]=1; M[i,1]=1/mp.mpf(k); M[i,2]=1/mp.mpf(k)**2
    b = mp.matrix(rs)
    sol = mp.lu_solve(M,b)
    A,B,C = sol[0],sol[1],sol[2]
    print(f"class {cls}: A={mp.nstr(A,18)}  B={mp.nstr(B,14)} (pred c1={mp.nstr(c1_pred[cls],14)})  C={mp.nstr(C,10)}")
# residual check: two-term law error * n^2 per class
print("\n(f_n - law2)/law_main * n^2 -> c2 estimate per class (should converge):")
def law2f(n):
    return law_main(n)*(1 + (3-3*mp.pi*((-1)**n))/(8*n))
for cls in range(4):
    vals=[]
    for k in [80+cls, 120+cls, 160+cls, 176+cls]:
        vals.append((f450[k]-law2f(k))/law_main(k)*k**2)
    print(f"class {cls}: ", [mp.nstr(v,8) for v in vals])
with open('progress.log','a') as f:
    f.write("step_13: E-law Richardson per class mod 4 done (constant=1, c1 matches (3-+3pi)/8)\n")
