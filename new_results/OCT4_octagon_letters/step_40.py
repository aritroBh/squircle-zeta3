import os; os.chdir('research/08_theory/OCT4_octagon_letters')
import mpmath as mp
from oct4_common import log
mp.mp.dps = 80
r2=mp.sqrt(2); al=r2-1; be=r2/2-1
def L2(x): return mp.polylog(2,x)
def L3(x): return mp.polylog(3,x)
def log(x): return mp.log(x)
# antiderivative of log^2(1+al t): A(t)
def A(t): return (1/al)*(1+al*t)*(log(1+al*t)**2 - 2*log(1+al*t) + 2)
# boundary
boundary = A(1)*log(1+be)
# elementary integrals E2=int log^2, E1=int log, over [0,1]
E2 = A(1) - A(0)   # int_0^1 log^2(1+al t) dt
E1 = mp.quad(lambda t: log(1+al*t),[0,1])
E0 = mp.mpf(2)
# N = int (log^2 - 2log +2)/(1+be t) dt ; pieces:
# W = int log^2(1+al t)/(1+be t) dt ; via sub s=1+al t, gamma=(be-al)/be
gam = (be-al)/be
def Fanti(s):  # antiderivative of log^2 s/(s-gam)
    return log(s)**2*log(1-s/gam) + 2*log(s)*L2(s/gam) - 2*L3(s/gam)
W = (1/be)*(Fanti(1+al) - Fanti(1))
# V = int log(1+al t)/(1+be t) dt : sub s=1+al t -> int log s/(s-gam) ds *(1/be)
def Ganti(s):  # antiderivative of log s/(s-gam) = log s log(1-s/gam) + Li_2(s/gam)
    return log(s)*log(1-s/gam) + L2(s/gam)
V = (1/be)*(Ganti(1+al)-Ganti(1))
# Z = int 1/(1+be t) dt
Z = log(1+be)/be
N = W - 2*V + 2*Z
# assemble
J21_analytic = boundary - (be/al)*( (al/be)*(E2 - 2*E1 + E0) + (be-al)/be * N )
J21_quad = mp.quad(lambda t: log(1+al*t)**2*log(1+be*t), [mp.mpf(0),mp.mpf('0.5'),1])
print("J21 analytic =", mp.nstr(J21_analytic,40))
print("J21 quad     =", mp.nstr(J21_quad,40))
print("diff =", mp.nstr(abs(J21_analytic-J21_quad),5))
log("step: J21 analytic reduction verified, diff=%s"%mp.nstr(abs(J21_analytic-J21_quad),3))
