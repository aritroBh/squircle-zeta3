
code = r'''
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, nstr, gamma, cos, sin
mp.dps = 60
def H(s): return sqrt(pi)/4*gamma(s/2)/gamma(s/2+mpf('1.5'))
rot = exp(mpc(0,pi/4))
def RePhi(x, N): return -pi*x/(2*sqrt(2)) + N*log(x) - N*log(1+sqrt(2)*x+x*x)/2
def RePhi2(x,N): return -N/x**2 + N*x*(x+sqrt(2))/(1+sqrt(2)*x+x*x)**2
def RePhi1(x,N): return -pi/(2*sqrt(2)) + N/x - N*(x+1/sqrt(2))/(1+sqrt(2)*x+x*x)
from mpmath import findroot

# For a range of N, compute: x0, RePhi(x0)+sqrt(pi N), |RePhi''(x0)|*sqrt(N), M3 := max|Phi'''| near, |omega0|*sqrt(N)
def Phi(x,N):
    s = rot*x
    return mpc(0,1)*pi*s/2 + N*log(s/(1+s))
def Phi1c(x,N):
    s = rot*x
    return mpc(0,1)*pi*rot/2 + N*rot/(s*(1+s))
def Phi2c(x,N):
    s = rot*x
    return -mpc(0,1)*N*(1+2*s)/(s*s*(1+s)**2)*rot**2
def Phi3c(x,N):
    s = rot*x
    return mpc(0,1)*N*(2+6*s+6*s**2)/(s**3*(1+s)**3)*rot**3
rows = []
for k in [50,100,150,200,250,300,350,400]:
    N = mpf(k+1)
    xs = sqrt(2*N/pi)
    x0 = findroot(lambda x: RePhi1(x,N), xs)
    R0 = RePhi(x0,N)
    w0 = abs(Phi1c(x0,N).imag)  # |omega0|
    p2 = abs(Phi2c(x0,N))
    re2 = abs(RePhi2(x0,N))
    M3 = abs(Phi3c(x0,N))
    rows.append((k,N,x0,R0,w0,p2,re2,M3))
    print("k=%3d: x0=%.5f (RePhi(x0)+th)*sqrtN=%.6f  w0*sqrtN=%.5f  |phi2|*sqrtN=%.5f  |RePhi2|*sqrtN=%.5f  M3*sqrtN=%.6f" % (
        k, x0, (R0+sqrt(pi*N))*sqrt(N), w0*sqrt(N), p2*sqrt(N), re2*sqrt(N), M3*sqrt(N)))

# check: RePhi(x) <= RePhi(x0) - c (x-x0)^2 on |x-x0| <= x0/2: estimate c = 1/2 min|RePhi2| on interval
print("\nconcavity curvature bound c(N) = 1/2 min|RePhi''| on [x0/2, 3x0/2]:")
for k in [100, 400]:
    N = mpf(k+1); xs = sqrt(2*N/pi)
    x0 = findroot(lambda x: RePhi1(x,N), xs)
    m = min(abs(RePhi2(x0/2 + (x0)*(mpf(j)/100), N)) for j in range(101))
    print("  k=%d: c=%s  c*sqrt(N)=%s" % (k, nstr(m/2,6), nstr(m/2*sqrt(N),6)))
# prefactor psi = rot*H(s)/s : |psi(x0)| * N^{5/4}
print("\nprefactor |psi(x0)| N^{5/4} and |psi'/psi| bounds:")
for k in [100, 400]:
    N = mpf(k+1); xs = sqrt(2*N/pi)
    x0 = findroot(lambda x: RePhi1(x,N), xs)
    s0 = rot*x0
    psi0 = rot*H(s0)/s0
    # |psi(x0+dl)/psi(x0)-1| <= Mpsi dl: estimate |psi'/psi| at x0
    h = mpf('1e-6')
    psip = (rot*H(rot*(x0+h))/(rot*(x0+h)) - psi0)/h
    print("  k=%d: |psi0| N^(5/4)=%s  |psi'/psi|=%s" % (k, nstr(abs(psi0)*N**mpf('1.25'),8), nstr(abs(psip/psi0),6)))
'''
open("step_14.py","w").write(code)
import subprocess, sys
r = subprocess.run([sys.executable,"step_14.py"],capture_output=True,text=True)
print(r.stdout); print("STDERR:", r.stderr[-1200:] if r.stderr else "")
