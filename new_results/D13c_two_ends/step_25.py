
code = r'''
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, nstr, cos, sin, gamma
import json
mp.dps = 80
R = {int(k): mpf(v) for k,v in json.load(open("mkmain_cache.json")).items()}
Mk = {int(k): mpf(v) for k,v in json.load(open("mkmain_cache.json")).items()}  # placeholder
def load_mk(fn):
    out = {}
    for line in open(fn):
        p = line.split("\t")
        if len(p) < 2: continue
        out[int(p[0])] = mpf(p[1].strip().split(" ")[0].strip("[]"))
    return out
Mkt = load_mk("../D13b_Mk_constant/mk_deg9.tsv")

# ---- q=2 removal: fit 4-param on [200,400], subtract, look for q=3 ----
def fitN(kmin,kmax,npar,wfun):
    ks = list(range(kmin,kmax+1))
    M = mp.matrix(len(ks),2*npar); y = mp.matrix(len(ks),1)
    for i,k in enumerate(ks):
        N = mpf(k+1); w2 = wfun(N); sc = exp(-w2)/N
        for j in range(npar):
            pw = N**(mpf('-0.5')*j)
            M[i,2*j] = sc*cos(w2)*pw; M[i,2*j+1] = sc*sin(w2)*pw
        y[i] = R[k]
    G = M.T*M; rhs = M.T*y
    sol = mp.lu_solve(G,rhs)
    return sol
sol2 = fitN(200,400,3,lambda N: sqrt(2*pi*N))
def q2fit(k):
    N = mpf(k+1); w2 = sqrt(2*pi*N); sc = exp(-w2)/N
    v = mpf(0)
    for j in range(3):
        pw = N**(mpf('-0.5')*j)
        v += sc*pw*(sol2[2*j]*cos(w2)+sol2[2*j+1]*sin(w2))
    return v
R2 = {k: R[k]-q2fit(k) for k in range(200,401)}
# fit R2 with rate sqrt(3 pi)
def fitCD_gen(seq, kmin, kmax, a):
    ks = list(range(kmin,kmax+1))
    M = mp.matrix(len(ks),2); y = mp.matrix(len(ks),1)
    for i,k in enumerate(ks):
        N = mpf(k+1); w = a*sqrt(N); sc = exp(-w)/N
        M[i,0] = sc*cos(w); M[i,1] = sc*sin(w); y[i] = seq[k]
    G = M.T*M; rhs = M.T*y
    sol = mp.lu_solve(G,rhs); res = M*sol-y
    return sol, sqrt(sum(res[i]**2 for i in range(len(ks)))/len(ks))
print("q=3 fit on R2 (rate sqrt(3 pi)), amplitude prediction 2 sqrt(pi)(1+1/sqrt3)=5.59 or 2 sqrt(pi)(1-1/sqrt3)=1.50:")
for kmin,kmax in [(200,400),(250,400),(300,400)]:
    sol,rms = fitCD_gen(R2,kmin,kmax,sqrt(3*pi))
    amp = sqrt(sol[0]**2+sol[1]**2); ph = mp.atan2(sol[1],sol[0])
    print("  [%3d,%3d]: amp=%s phase=%s rms=%s" % (kmin,kmax,nstr(amp,8),nstr(ph,6),nstr(rms,2)))
print("  rate test for R2:")
for name,a in [("sqrt(3 pi)",sqrt(3*pi)),("2 sqrt(pi)",2*sqrt(pi)),("sqrt(2 pi)",sqrt(2*pi))]:
    sol,rms = fitCD_gen(R2,250,400,a)
    print("   %s: amp=%s rms=%s" % (name,nstr(sqrt(sol[0]**2+sol[1]**2),6),nstr(rms,2)))

# ---- kappa1 check: brief's kappa1 vs D13b beta_+ ----
z = mpc(0,pi/2)
beta_p = 1/(2*sqrt(z)) + sqrt(z)/2 - z**mpf('1.5')/12
k1_brief = ((24+12*pi+pi**2) + mpc(0,1)*(24-12*pi+pi**2))/(48*sqrt(pi))
print("\nkappa1 brief =", nstr(k1_brief, 15))
print("beta_+ (D13b) =", nstr(beta_p, 15))
print("conj(beta_+) =", nstr(beta_p.conjugate(), 15))
print("|brief - conj(beta_+)| =", nstr(abs(k1_brief-beta_p.conjugate()),4))
c_r = (1/(2*sqrt(pi)) + sqrt(pi)/4 + pi**mpf('1.5')/48)
c_i = (-1/(2*sqrt(pi)) + sqrt(pi)/4 - pi**mpf('1.5')/48)
print("c_r, c_i =", nstr(c_r,12), nstr(c_i,12))

# ---- D11 jump recomputation WITH 1/eta Laplace prefactor ----
# B_F(t) ~ c_+ (1-t/A)^{1/2}, c_+ = -sqrt(2 pi) e^{3 i pi/4}; Disc across cut: 2 i c_+ (t/A-1)^{1/2}
# J(eta) = (1/eta) int_cut e^{-t/eta} Disc B dt  ~  (1/eta) e^{-A/eta} 2 i c_+ int_0^inf e^{-x/eta} (x/A)^{1/2} dx
#        = (1/eta) e^{-A/eta} 2 i c_+ A^{-1/2} Gamma(3/2) eta^{3/2} = 2 i c_+ A^{-1/2} (sqrt pi/2) eta^{1/2} e^{-A/eta}
cp = -sqrt(2*pi)*exp(mpc(0,3*pi/4)); A = mpc(0,pi/2)
Jc = 2*mpc(0,1)*cp*A**mpf('-0.5')*sqrt(pi)/2
print("\ncorrect jump J =", nstr(Jc,8), "* eta^{1/2} e^{-A/eta}   (magnitude, up to cut-orientation sign)")
print("  |coef| =", nstr(abs(Jc),10), " vs 2 sqrt pi =", nstr(2*sqrt(pi),10))
print("  D11 jump: 2 sqrt(pi) e^{3i pi/4} eta^{3/2}: ratio D11/correct at eta=i eps:", 
      nstr(2*sqrt(pi)*exp(mpc(0,3*pi/4))*(mpc(0,1))**mpf('1.5')/(Jc*(mpc(0,1))**mpf('0.5')),8))
print("  = should be i*... i.e. exactly one power of eta off")
'''
open("step_11.py","w").write(code)
import subprocess, sys
r = subprocess.run([sys.executable,"step_11.py"],capture_output=True,text=True)
print(r.stdout); print("STDERR:", r.stderr[-1200:] if r.stderr else "")
