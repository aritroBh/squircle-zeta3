import subprocess, os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
code = r'''
import mpmath as mp
mp.mp.dps = 120
from math import factorial

N = 130
Zv = {n: mp.zeta(n) for n range(2, N+3)} if False else {n: mp.zeta(n) for n in range(2, N+3)}
def series_mul(a, b, N):
    c = [mp.mpf(0)]*(N+1)
    for i in range(N+1):
        if a[i]==0: continue
        for j in range(N+1-i):
            if b[j]!=0: c[i+j]+=a[i]*b[j]
    return c
def series_exp(a, N):
    out = [mp.mpf(1)]+[mp.mpf(0)]*N
    for m in range(1, N+1):
        out[m] = sum(j*a[j]*out[m-j] for j in range(1, m+1) if a[j]!=0)/m
    return out
def rhombus_fn(q, N):
    a = [(1 if n%2==0 else q)/mp.factorial(n) for n in range(N+1)]
    g = [mp.mpf(0)]*(N+1); g[0]=mp.mpf(1)
    for k in range(1, N+1):
        g[k] = (a[k] - sum(g[i]*g[k-i] for i in range(1,k)))/2
    mu = {}
    mu[0] = [mp.mpf(0)]*(N+1)
    for n in range(1, N+1): mu[0][n] = ((-1)**n)/n
    for n in range(2, N+1):
        mu[0][n] += ((-1)**n)*Zv[n]/n*(mp.mpf(2)**(1-n)-1)
    for j in range(2, N+1, 2):
        arr = [mp.mpf(0)]*(N+1)
        for n in range(max(2,j), N+1):
            arr[n-j] = 2*((-1)**n)*Zv[n]*factorial(n-1)/factorial(n-j)/mp.mpf(2)**(n-j)
        mu[j] = arr
    A0 = series_exp(mu[0], N)
    e = [[mp.mpf(1)]+[mp.mpf(0)]*N]
    for k in range(1, N+1):
        ek = [mp.mpf(0)]*(N+1)
        for j in range(2, k+1, 2):
            conv = series_mul(mu[j], e[k-j], N)
            for m in range(N+1):
                ek[m] += conv[m]/factorial(j-1)
        e.append([v/k for v in ek])
    dsA = {k: [factorial(k)*v for v in series_mul(A0, e[k], N)] for k in range(0, N+1)}
    I = [mp.mpf(0)]*(N+1)
    for k in range(0, N+1, 2):
        gk = g[k]
        if gk==0: continue
        ser = dsA[k]
        for m in range(N+1-k):
            I[k+m] += gk*ser[m]
    return I

def test2(A, B, nmax=125):
    q = mp.mpf(A*A-B*B)/mp.mpf(A*A+B*B)
    r = mp.sqrt(A*A+B*B); L = mp.log(mp.mpf(B)/A); rho = mp.sqrt(L**2+mp.pi**2/4)
    th0 = mp.atan2(mp.pi/2, L)
    f = rhombus_fn(q, N)
    zlist = [L + 1j*mp.pi/2, -L + 1j*mp.pi/2, L - 1j*mp.pi/2, -L - 1j*mp.pi/2]
    clist = [A*mp.sqrt(2)/r, B*mp.sqrt(2)/r, A*mp.sqrt(2)/r, B*mp.sqrt(2)/r]
    def bracket(n):
        return B*((-1)**n)*mp.cos((n-mp.mpf('0.5'))*th0) - A*mp.sin((n-mp.mpf('0.5'))*th0)
    def law2(n):
        tot = mp.mpc(0)
        for z, c in zip(zlist, clist):
            tot += c*mp.sqrt(-z)*z**(-n)*(1 + mp.mpf(3)/(8*n) - 3*z/(4*n))
        return (-1/(2*mp.sqrt(mp.pi)))*n**(-mp.mpf('1.5'))*tot
    errs = []
    for n in range(60, nmax+1):
        if abs(bracket(n)) < 0.6: continue
        lv = law2(n)
        errs.append((n, abs(1 - (f[n]/mp.factorial(n))/lv.real)))
    w = max(e for _, e in errs)
    print(f"A={A},B={B}: two-term law, safe n={len(errs)}, worst |1-f_n/(n! law2H)| = {mp.nstr(w,4)}")
    print(f"   at n={errs[-1][0]}: {mp.nstr(errs[-1][1],4)}  x n^2 = {mp.nstr(errs[-1][1]*errs[-1][0]**2,5)}")
    # also 1-term law worst
    def law1(n):
        return -(mp.sqrt(2)/(r*mp.sqrt(mp.pi)))*rho**(mp.mpf('0.5')-n)*n**(-mp.mpf('1.5'))*bracket(n)
    errs1 = []
    for n in range(60, nmax+1):
        if abs(bracket(n)) < 0.6: continue
        errs1.append(abs(1-(f[n]/mp.factorial(n))/law1(n)))
    print(f"   one-term law worst |1-ratio| = {mp.nstr(max(errs1),4)}")

test2(3, 2)
test2(1, 3)
'''
with open(wd+"/step_43.py","w") as fp: fp.write(code)
t0=time.time()
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_43.py"], capture_output=True, text=True, cwd=wd, timeout=1700)
print(r.stdout[-2500:]); print("STDERR:", r.stderr[-800:]); print("elapsed", time.time()-t0)
