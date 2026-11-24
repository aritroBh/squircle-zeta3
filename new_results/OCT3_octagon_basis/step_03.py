import json, time
from mpmath import mp, mpf, sqrt, sin, cos, pi, log, exp, quad, nstr, factorial

mp.dps = 210
t0=time.time()
N = 4  # Taylor order in t = p-1

def s_mul(a,b):
    c=[mpf(0)]*N
    for i in range(N):
        ai=a[i]
        for j in range(N-i):
            c[i+j]+=ai*b[j]
    return c
def s_inv(a):
    c=[mpf(0)]*N
    ia=1/a[0]; c[0]=ia
    for k in range(1,N):
        s=mpf(0)
        for i in range(1,k+1): s+=a[i]*c[k-i]
        c[k]=-s*ia
    return c
def s_div(a,b): return s_mul(a, s_inv(b))
def s_exp(a):
    c=[mpf(0)]*N; c[0]=exp(a[0])
    for k in range(1,N):
        s=mpf(0)
        for i in range(1,k+1): s+=i*a[i]*c[k-i]
        c[k]=s/k
    return c
def s_log(a):
    c=[mpf(0)]*N; c[0]=log(a[0])
    ia=s_inv(a)
    for k in range(1,N):
        s=mpf(0)
        for i in range(0,k): s+=(i+1)*a[i+1]*ia[k-1-i]
        c[k]=s/k
    return c
def s_sqrt(a):
    h=mpf('0.5'); return s_exp([h*x for x in s_log(a)])

ONE_T=[mpf(1),mpf(1),mpf(0),mpf(0)]
PI4=pi/4

def fseries(th):
    f1=cos(th-PI4); d1=-sin(th-PI4)
    f2=sin(th);     d2=cos(th)
    f3=cos(th+PI4); d3=-sin(th+PI4)
    f4=cos(th);     d4=-sin(th)
    u=[mpf(0)]*N; w=[mpf(0)]*N
    for f,d in ((f1,d1),(f2,d2),(f3,d3),(f4,d4)):
        L=log(f)
        ek=(mpf(1), L, L*L/2, L*L*L/6)
        for k in range(N):
            u[k]+=f*ek[k]; w[k]+=d*ek[k]
    s   = [-x for x in s_div(s_log(u), ONE_T)]
    sth = [-x for x in s_div(w, u)]
    r   = s_exp(s)
    rth = s_mul(r, sth)
    rr  = s_mul(r,r); rt  = s_mul(rth,rth)
    return s_sqrt([rr[i]+rt[i] for i in range(N)])

Pk=[]
for k in range(4):
    val = 8*factorial(k)*quad(lambda th: fseries(th)[k], [0, PI4])
    Pk.append(val)
    print(f"P^{(k)}(1) computed:", nstr(val, 40), f"  [{time.time()-t0:.1f}s]")

data = json.load(open('research/08_theory/D15b_polygon_law/vals_m2_8_200.json'))
file4 = [mpf(x) for x in data["4"]]
print("\n--- comparison with file key '4' ---")
for k in range(4):
    diff = Pk[k]-file4[k]
    print(f"k={k}: |computed-file| = {nstr(abs(diff),5)}")
print("\nexact P(1) check:", nstr(Pk[0]-8*(sqrt(2)-1)*sqrt(2-sqrt(2)), 5))
