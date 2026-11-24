code = '''"""RST1R step_11: R11 / Item 1 (numerical by design) independent checks.
(a) sympy: A_k^[2] identity  -s^k/(s+2)^{k+1} - s^{k-1}/(s+2)^k = -2 s^{k-1}(s+1)/(s+2)^{k+1};
    saddle equation (i pi/2) s(s+1)(s+2) = s^2 - 2(k-1)s - 2(k-1)  from Phi'(s)=0.
(b) rest/R_k at k=100,150,...,400 (exact M_k balls vs mpmath quadratures).
(c) saddle constants |g_k| N e^{th2} -> sqrt pi, arg-th2 -> pi/2 at k=400,2000,1e4,5e4;
    Re Phi(s*) expansion: Re Phi + sqrt(2 pi (k-1)) + (1/2) log(4(k-1)/pi) -> 0.
(d) subdominance |R_q3|/|R_q2|, |R_KE3|/|R_q2| at k=100,200,400.
"""
import mpmath as mp, re, sympy as sp
mp.mp.dps = 80

s, k = sp.symbols('s k')
expr1 = -s**k/(s+2)**(k+1) - s**(k-1)/(s+2)**k + 2*s**(k-1)*(s+1)/(s+2)**(k+1)
print("(a1) A_k^[2] identity difference:", sp.simplify(expr1))
assert sp.simplify(expr1) == 0
# Phi'(s) = i pi/2 + (k-1)/s + 1/(s+1) - (k+1)/(s+2); multiply by s(s+1)(s+2)
Ph = sp.I*sp.pi/2 + (k-1)/s + 1/(s+1) - (k+1)/(s+2)
poly = sp.expand(Ph * s*(s+1)*(s+2))
print("(a2) Phi'*s(s+1)(s+2) =", poly)
claimed = sp.I*sp.pi/2*s*(s+1)*(s+2) - (s**2 - 2*(k-1)*s - 2*(k-1))
print("(a2) difference from claimed:", sp.expand(poly - claimed))
assert sp.expand(poly - claimed) == 0
print("(a) algebraic identities: PROVED")

pi=mp.pi; j=1j; omega=mp.expjpi(mp.mpf(1)/4)
def H(x): return mp.beta(x/2, mp.mpf(3)/2)/2
sq2pi=mp.sqrt(2*pi); sqpi=mp.sqrt(pi)
rows={}
for line in open("<project>/research/08_theory/D13f_sweep/mk_deg9.tsv"):
    p=line.rstrip("\\n").split("\\t")
    if p[0].isdigit(): rows[int(p[0])]=mp.mpf(re.match(r"\\[([-0-9.eE+]+)",p[1]).group(1))
def Mmain(k):
    N=k+1; rs=mp.sqrt(2*N/pi)
    f=lambda r: mp.expj(pi*(r*omega)/2)*(r*omega/(1+r*omega))**N*H(r*omega)/(r*omega)*omega
    return 2/pi*mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]).imag
def Rmodel2(k):
    rs=2*mp.sqrt(k/pi)
    def f(r):
        s=r*omega
        return H(s)*mp.expj(pi*s/2)*s**(k-1)*(s+1)/(s+2)**(k+1)*omega
    return -4/pi*mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]).imag

print("(b) rest/R_k (proof: 0.99971,1.00262,1.00004,0.99981,0.99999,0.99997,0.99999)")
for kk in [100,150,200,250,300,350,400]:
    rest=rows[kk]-Mmain(kk); Rm=Rmodel2(kk)
    print(f"  k={kk}: rest/R_k = {mp.nstr(rest/Rm,7)}", flush=True)

def saddle(kk):
    r0=2*mp.sqrt((kk-1)/pi)
    dPhi=lambda x: j*pi/2+(kk-1)/x+1/(x+1)-(kk+1)/(x+2)
    Phi=lambda x: j*pi*x/2+(kk-1)*mp.log(x)+mp.log(x+1)-(kk+1)*mp.log(x+2)
    ddPhi=lambda x: -(kk-1)/x**2-1/(x+1)**2+(kk+1)/(x+2)**2
    s_=mp.findroot(dPhi, r0*omega)
    return s_, Phi, ddPhi
print("(c) saddle constants (proof: 1.98950,1.86506,1.81302,1.79043 -> sqrt pi; args 1.51096,1.54405,1.55884,1.56545 -> pi/2)")
for kk in [400,2000,10000,50000]:
    N=kk+1; th=sq2pi*mp.sqrt(N)
    s_,Phi,ddPhi = saddle(kk)
    g=-(4/pi)*omega*mp.exp(Phi(s_))*H(s_)*mp.sqrt(2*pi/(-ddPhi(s_)))
    resid = mp.re(Phi(s_)) + sq2pi*mp.sqrt(kk-1) + mp.mpf('0.5')*mp.log(4*(kk-1)/pi)
    print(f"  k={kk}: |g|N e^th2 = {mp.nstr(abs(g)*N*mp.e**th,8)}; arg-th2 = {mp.nstr(mp.fmod(mp.arg(g)-th+pi,2*pi)-pi,8)}; Re-Phi resid = {mp.nstr(resid,6)}", flush=True)

def Rmodel3(k):
    rs=mp.sqrt(6*k/pi)
    def f(r):
        s=r*omega
        return H(s)*mp.expj(pi*s/2)*s**(k-2)/(s+3)**(k-1)*omega
    return (1/pi)*mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]).imag
def RKE3(k):
    rs=mp.sqrt(2*k/(3*pi))
    def f(r):
        s=r*omega
        return H(s)*mp.expj(3*pi*s/2)*s**k/(1+s)**(k+1)*omega
    return (2/pi)*mp.quad(f,[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]).imag
print("(d) subdominance (proof: q3 9.8e-4,4.8e-5,4.0e-6; KE 4.5e-3,2.5e-4,3.4e-5)")
for kk in [100,200,400]:
    r2=abs(Rmodel2(kk)); r3=abs(Rmodel3(kk)); rke=abs(RKE3(kk))
    print(f"  k={kk}: |R_q3|/|R_q2| = {mp.nstr(r3/r2,6)}, |R_KE3|/|R_q2| = {mp.nstr(rke/r2,6)}", flush=True)
print("R11 numerical checks complete.")
'''
open('step_11.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_11.py'],capture_output=True,text=True,timeout=3000)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
