# Compare M_k from two Cauchy circles; growth |M_k|^(1/k); compare with leading law of eq:Mk.
import mpmath as mp
mp.mp.dps = 40
def load(f): return {int(l.split()[0]): mp.mpf(l.split()[1]) for l in open(f) if l.strip()}
A = load("cauchy_r095.log"); B = load("cauchy_r090.log"); C = load("cauchy_r025.log")
s2=mp.sqrt(2); a=1+s2; b=s2-1; L=mp.log(a); lam=mp.log(2)
M1c = (19*s2/48-mp.mpf(1)/24)*mp.pi**2+lam**2/2+2*s2*lam*L-(1+5*s2/2)*L**2-5*s2*mp.re(mp.polylog(2,b))
mp.mp.dps=110
s2=mp.sqrt(2); a=1+s2; b=s2-1; L=mp.log(a); lam=mp.log(2)
M1c = (19*s2/48-mp.mpf(1)/24)*mp.pi**2+lam**2/2+2*s2*lam*L-(1+5*s2/2)*L**2-5*s2*mp.re(mp.polylog(2,b))
M1def = mp.mpf(open("cauchy_r025.log").read().split("\n")[1].split()[1])
d = abs(M1def-M1c); print("M1 (definition, Cauchy r=1/4, N=160, dps100) vs closed form: diff =", mp.nstr(d,5), " -> agreeing digits ~", int(-mp.log10(d/abs(M1c))))
M0c = lam - s2*L; M0def = mp.mpf(open("cauchy_r025.log").read().split("\n")[0].split()[1]); print("M0 diff:", mp.nstr(abs(M0def-M0c),5))
mp.mp.dps=40
print(" k   M_k(r=.95)         rel.diff(r=.9)   |M_k|^(1/k)   law -(sqrt(pi)/N)e^{-sqrt(pi N)}cos(sqrt(pi N))   ratio")
for k in [10,20,40,60,80,100,120,140,160,180,200,220]:
    N=k+1; law = -(mp.sqrt(mp.pi)/N)*mp.exp(-mp.sqrt(mp.pi*N))*mp.cos(mp.sqrt(mp.pi*N))
    print(k, mp.nstr(A[k],12), mp.nstr(abs(A[k]-B[k])/abs(A[k]),3), mp.nstr(abs(A[k])**(mp.mpf(1)/k),8), mp.nstr(law,8), mp.nstr(A[k]/law,6))
# running max of |M_k|^(1/k) over windows and comparison with exp(-sqrt(pi/k))
print("window max of |M_k|^(1/k) vs exp(-sqrt(pi/k)) (envelope prediction):")
for lo in range(20,221,40):
    w = [abs(A[k])**(mp.mpf(1)/k) for k in range(lo, min(lo+40,221))]
    print(lo, "-", min(lo+39,220), mp.nstr(max(w),8), mp.nstr(mp.exp(-mp.sqrt(mp.pi/(lo+20))),8))
