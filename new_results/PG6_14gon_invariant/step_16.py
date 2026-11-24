import mpmath as mp, pickle, sys
sys.path.insert(0,'.')
import ftfield as F
mp.mp.dps = 80
m = 7
def fsum(th, p):
    return sum(mp.exp(p*mp.log(abs(mp.cos(th - i*mp.pi/m)))) for i in range(m))
def r_and_rp(th, p):
    f = fsum(th,p)
    fp = sum(-mp.sin(th - i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    return f**(-1/p), -(1/p)*f**(-1/p-1)*fp
def integrand(th, p):
    r, rp = r_and_rp(th,p); return mp.sqrt(r*r+rp*rp)
def perim_chamber(p, lo, hi):
    return 2*m*mp.quad(lambda th: integrand(th,p), [lo,hi], method='tanh-sinh', maxdegree=8)
c3 = { -4: -mp.mpf(7)/240, -3: mp.mpf(3)/10, -2: -mp.mpf(169)/120, -1: mp.mpf(61)/30,
       1: -mp.mpf(61)/30, 2: mp.mpf(169)/120, 3: -mp.mpf(3)/10, 4: mp.mpf(7)/240 }
def P3_fd(lo,hi,h):
    vals = {k: perim_chamber(1+k*h,lo,hi) for k in range(-4,5)}
    return sum(c3[k]*vals[k] for k in c3)/h**3
h = mp.mpf('0.0025')
# (a) gold chamber WITH the 1e-9 slivers removed -> should reproduce gold value
lo_g = mp.pi/2+mp.mpf('1e-9'); hi_g = mp.pi/2+mp.pi/m-mp.mpf('1e-9')
print("gold chamber (slivers removed) P3 =", mp.nstr(P3_fd(lo_g,hi_g,h),30))
# (b) FULL gold chamber [pi/2, pi/2+pi/7] (no sliver removal) -> should give truth
lo_f = mp.pi/2; hi_f = mp.pi/2+mp.pi/m
print("full offset chamber [pi/2,pi/2+pi/7] P3 =", mp.nstr(P3_fd(lo_f,hi_f,h),30))
# (c) centered chamber (truth)
print("centered chamber P3 =", mp.nstr(P3_fd(-mp.pi/14,mp.pi/14,h),30))
print("truth (chain) = -5.5080836980534783")
print("gold          = -5.5080732424212877")
open('progress.log','a').write('PG6 s_A3: gold chamber(sliver) reproduces gold; full/centered give truth. Bug = 1e-9 vertex slivers.\n')
