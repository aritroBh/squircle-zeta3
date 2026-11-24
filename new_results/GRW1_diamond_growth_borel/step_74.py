import subprocess, os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
code = r'''
import mpmath as mp
mp.mp.dps = 120
from math import factorial

N = 90
Zv = {n: mp.zeta(n) for n in range(2, N+3)}
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
    for n in range(1, N+1):
        mu[0][n] = ((-1)**n)/n
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
    return I, g

def test_case(A, B, Nmax=85):
    q = mp.mpf(A*A-B*B)/mp.mpf(A*A+B*B)
    r = mp.sqrt(A*A+B*B); L = mp.log(mp.mpf(B)/A); rho = mp.sqrt(L**2+mp.pi**2/4)
    I, g = rhombus_fn(q, N)
    th0 = mp.atan2(mp.pi/2, L)
    def bracket(n):
        return B*((-1)**n)*mp.cos((n-mp.mpf('0.5'))*th0) - A*mp.sin((n-mp.mpf('0.5'))*th0)
    def law(n):
        return -(mp.sqrt(2)/(r*mp.sqrt(mp.pi)))*rho**(mp.mpf('0.5')-n)*n**(-mp.mpf('1.5'))*bracket(n)*mp.factorial(n)
    f = I
    safe = [(n, f[n]/law(n)) for n in range(40, Nmax+1) if abs(bracket(n)) > 0.6]
    wr = max(abs(v-1) for _, v in safe)
    print(f"A={A},B={B}: q={mp.nstr(q,6)} rho={mp.nstr(rho,12)} th0={mp.nstr(th0,10)}")
    print(f"   f_1={mp.nstr(f[1],6)} (=-1), safe n count {len(safe)}, worst |ratio-1| = {mp.nstr(wr,4)}")
    for n, v in safe[-4:]:
        print(f"   n={n}: ratio={mp.nstr(v, 16)}")
    bn = abs(f[Nmax]/mp.factorial(Nmax))**(mp.mpf(1)/Nmax)
    print(f"   |b_{Nmax}|^(1/n) = {mp.nstr(bn, 12)}   1/rho = {mp.nstr(1/rho, 12)}")

test_case(3, 2)   # q = 5/13
test_case(1, 3)   # q = -4/5
'''
with open(wd+"/step_42.py","w") as fp: fp.write(code)
t0=time.time()
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_42.py"], capture_output=True, text=True, cwd=wd, timeout=1700)
print(r.stdout[-3000:]); print("STDERR:", r.stderr[-800:]); print("elapsed", time.time()-t0)
