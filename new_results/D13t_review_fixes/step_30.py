with open('progress.log','a') as f:
    f.write("turn10-11 R3 (H1): saddle branch rho(t)=sqrt(1+t^2/(1+rt2 t)) solves Re Phi'(1/t)=0 (resid 9.9e-48, 300 pts);\n"
            "  branch eq rho^2(1+rt2 t)=1+rt2 t+t^2 resid 3.2e-58. B_true exact form: worst rel resid 1.01e-32 over 5000 pts t in [1e-14,1e3].\n"
            "  phi = a(1-2rho)/t + N atan(zz) resid 5.5e-47. zz<1 for ALL t (max 0.9986 at t=1e3) => alternating-series bracket valid everywhere;\n"
            "  A1n bracket never fails: min(A1n-L)=1.79e-50 (100dps; the 60dps -6.4e-38 was noise), min(U-A1n)=5e-26, t in [1e-12,1e8], 20001 pts.\n"
            "  arb Btrue_box well-formed, contains truth at 100000/100000 pts (max viol -1.48e-6, matches D13s 1.5e-6); B_abs=0.16666814<=breg_abs=0.75.\n"
            "  R3 VERDICT: HOLDS (exact form verified; bracket proved via |zz|<1 all t; no failure found).\n")
import mpmath as mp
mp.mp.dps = 50
mpr2=mp.sqrt(2); mpp=mp.pi; sqpi=mp.sqrt(mpp); mbb=mp.e**(1j*mpp/4); ma=mpp/(2*mpr2)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(mbb*x)/x
# (i) x^2|g(x)| <= 1 on (0,1], dense, incl tiny x
mx1 = mp.mpf('0'); at=None
for k in range(1,200001):
    x = mp.mpf(k)/200000
    v = x*x*abs(g(x))
    if v > mx1: mx1=v; at=x
mxtiny = max((mp.mpf(10)**(-k))**2*abs(g(mp.mpf(10)**(-k))) for k in range(1,31))
print(f"R4 (i) max x^2|g(x)| on (0,1] (200000 pts) = {mp.nstr(mx1,12)} at x={float(at)}; tiny-x max = {mp.nstr(mxtiny,12)} (<=1)")
# series identity + c_k properties
def c_k(k): return (-1)**(k+1)*mp.binomial(mp.mpf('0.5'), k)
print(f"    c_k>0 for k=1..60: {all(c_k(k)>0 for k in range(1,61))}; sum c_k k=1..2e4 = {mp.nstr(mp.nsum(c_k,[1,20000]),10)} (=1?)")
resid = mp.mpf('0')
for xr in [mp.mpf('0.001'), mp.mpf('0.05'), mp.mpf('0.37'), mp.mpf('0.83'), mp.mpf('1.0')]:
    s = xr*mbb
    r = s*Hc(s) - (1 - mp.nsum(lambda kk: c_k(int(kk))*s/(s+2*kk), [1, mp.inf], tol=mp.mpf('1e-40')))
    resid = max(resid, abs(r))
print(f"    series identity sHc(s)=1-sum c_k s/(s+2k): worst resid at 5 ray pts = {mp.nstr(resid,6)}")
# beta-integral identity Hc(s) = (1/2) int_0^1 u^{s/2-1}(1-u)^{1/2} du
s = mp.mpf('0.37')*mbb
be = mp.mpf('0.5')*mp.quad(lambda u: u**(s/2-1)*(1-u)**mp.mpf('0.5'), [0,1])
print(f"    beta-integral identity resid = {mp.nstr(abs(Hc(s)-be),6)}")
