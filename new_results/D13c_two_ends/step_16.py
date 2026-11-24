
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, nstr, gamma, cos, sin
import json
mp.dps = 70
def load_mk(fn):
    out = {}
    for line in open(fn):
        p = line.split("\t")
        if len(p) < 2: continue
        out[int(p[0])] = mpf(p[1].strip().split(" ")[0].strip("[]"))
    return out
Mkt = load_mk("../D13b_Mk_constant/mk_deg9.tsv")
def H(s): return sqrt(pi)/4*gamma(s/2)/gamma(s/2+mpf("1.5"))
def sstar(N): return -mpf("0.5") + sqrt(mpf("0.25") + mpc(0,2*N/pi))
def phi_s(s,N): return mpc(0,1)*pi*s/2 + N*log(s/(1+s))
def phi2_s(s,N): return -mpc(0,1)*N*(1+2*s)/(s*s*(1+s)**2)
def SV(k):
    N = mpf(k+1)
    s = sstar(N)
    return exp(phi_s(s,N))*sqrt(2*pi/(-phi2_s(s,N)))*H(s)/s
def Mkmain_of(k):  # from cache
    R = {int(kk): mpf(v) for kk,v in json.load(open("mkmain_cache.json")).items()}
    return Mkt[k]-R[k]
print("k: |SV - lead|/(N^-3/2 e^-th), |Mmain - SV|/(N^-3/2 e^-th), |Mmain-lead|/(N^-3/2 e^-th)")
for k in [50,100,150,200,250,300,350,400]:
    N = mpf(k+1); th = sqrt(pi*N)
    sv = SV(k); mm = Mkmain_of(k)
    lead = -(sqrt(pi)/N)*exp(-th)*cos(th)
    sc = N**mpf("-1.5")*exp(-th)
    e1 = abs((2/pi)*sv.imag - lead)/sc          # saddle value vs leading law
    e2 = abs(mm - (2/pi)*sv.imag)/sc            # Laplace remainder vs saddle value
    e3 = abs(mm - lead)/sc                      # total
    print("  k=%3d:  SV-lead: %s   Mmain-SV: %s   total: %s" % (k, nstr(e1,5), nstr(e2,5), nstr(e3,5)))
# fit kappa1' for M_k^main:  Mmain = -(sqrt pi/N)e^{-th} Re[e^{i th}(1+K1 N^-1/2+K2/N+K3 N^-3/2)]
print("\nfit K1,K2,K3 of M_k^main on [100,400] and compare with brief kappa:")
def fitK(kmin,kmax):
    ks = list(range(kmin,kmax+1))
    M = mp.matrix(len(ks),6); y = mp.matrix(len(ks),1)
    for i,k in enumerate(ks):
        N = mpf(k+1); th = sqrt(pi*N)
        Yk = -Mkmain_of(k)*N*exp(th)/sqrt(pi)
        base = mpc(cos(th), sin(th))
        y[i] = Yk
        for j in range(3):
            pw = N**(mpf("-0.5")*(j+1))
            M[i,2*j]   = cos(th)*pw; M[i,2*j+1] = -sin(th)*pw
            if j==0 and False: pass
        M[i,0] = cos(th)*N**mpf("-0.5"); M[i,1] = -sin(th)*N**mpf("-0.5")
        # we want Yk = cos + Re[e^{i th}(K1 N^-1/2+...)]: unknowns Re Kj, Im Kj
    # unknowns: [ReK1, ImK1, ReK2, ImK2, ReK3, ImK3]; Yk - cos(th) = Re[e^{i th} sum] = cos*ReK - sin*ImK times pw
    for i,k in enumerate(ks):
        N = mpf(k+1); th = sqrt(pi*N)
        y[i] = -Mkmain_of(k)*N*exp(th)/sqrt(pi) - cos(th)
        for j in range(3):
            pw = N**(mpf("-0.5")*(j+1))
            M[i,2*j] = cos(th)*pw; M[i,2*j+1] = -sin(th)*pw
    G = M.T*M; rhs = M.T*y
    sol = mp.lu_solve(G,rhs); res = M*sol-y
    return sol, sqrt(sum(res[i]**2 for i in range(len(ks)))/len(ks))
k1b = ((24+12*pi+pi**2) + mpc(0,1)*(24-12*pi+pi**2))/(48*sqrt(pi))
print("brief kappa1 =", nstr(k1b,10), " kappa2 = 0.852808-0.215617i kappa3=0.948724-0.426650i")
for kmin,kmax in [(100,400),(200,400),(300,400)]:
    sol,rms = fitK(kmin,kmax)
    K1 = sol[0]+mpc(0,1)*sol[1]; K2 = sol[2]+mpc(0,1)*sol[3]; K3 = sol[4]+mpc(0,1)*sol[5]
    print("  [%3d,%3d]: K1=%s (err %s) K2=%s K3=%s rms=%s" % (kmin,kmax,nstr(K1,8),nstr(abs(K1-k1b),3),nstr(K2,7),nstr(K3,7),nstr(rms,2)))
