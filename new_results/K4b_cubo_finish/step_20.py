import os, json
os.chdir("research/08_theory/K4b_cubo_finish")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
EPS = mpf('0.005')
out = json.load(open("fd_Svals.json"))
S = {int(k): mpf(v) for k, v in out.items()}
Sdk = [mpf(s) for s in json.load(open("tetra_Sdk_final.json"))["Sdk"]]

def fornberg_weights(N, k, eps):
    # weights w_j, j=-N..N, exact for degree <= 2N: sum w_j f(j eps) = f^{(k)}(0)
    idx = list(range(-N, N+1))
    A = mp.matrix(2*N+1, 2*N+1)
    for r in range(2*N+1):
        for c in range(2*N+1):
            A[r, c] = (idx[c]*eps)**r
    b = mp.matrix(2*N+1, 1)
    b[k] = mp.factorial(k)
    w = mp.lu_solve(A, b)
    return {idx[c]: w[c] for c in range(2*N+1)}

def deriv(N, k):
    w = fornberg_weights(N, k, EPS)
    return sum(wj*S[j] for j, wj in w.items()), sum(abs(wj) for j, wj in w.items())

print(f"{'k':>2} {'stencil':>8} {'FD derivative':>42} {'err vs series·k!':>18} {'Σ|w|·2.5e-32':>14}")
results = {}
for k in [0, 1, 2, 3, 4]:
    target = mp.factorial(k)*Sdk[k]
    for N in [6, 5, 4, 3]:
        v, sw = deriv(N, k)
        err = abs(v - target)
        noise = sw*mpf('2.5e-32')
        if N == 6: results[k] = v
        print(f"{k:>2} {2*N+1:>8} {mp.nstr(v,30):>42} {mp.nstr(err,2):>18} {mp.nstr(noise,2):>14}")
    print(f"   target k!·[d^{k}]S = {mp.nstr(target,40)}")

print("\n=== SUMMARY: FD (13-pt) vs series ===")
fd_res = {}
for k in range(5):
    target = mp.factorial(k)*Sdk[k]
    err = abs(results[k]-target)
    fd_res[f"d{k}"] = {"FD": mp.nstr(results[k],50), "series_k!": mp.nstr(target,50), "err": mp.nstr(err,3)}
    agree = -mp.log10(err/abs(target)) if err>0 else mpf(99)
    print(f"k={k}: agreement {mp.nstr(agree,4)} digits   [{ 'PASS' if agree>=15 else 'CHECK'}]")
json.dump(fd_res, open("fd_check.json","w"), indent=1)
print("saved fd_check.json")
