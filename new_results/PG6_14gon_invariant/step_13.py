import mpmath as mp, pickle, sys
sys.path.insert(0,'.')
import ftfield as F
mp.mp.dps = 80
m = 7
g = pickle.load(open('geom14.pkl','rb'))
LL = mp.sqrt(F.feval(g['L'],80))   # edge length L
hh = F.feval(g['h'],80)            # apothem h
print("edge length L =", mp.nstr(LL,40))
print("apothem h     =", mp.nstr(hh,40))
print("14*L          =", mp.nstr(14*LL,40))

def fsum(th, p):
    return sum(mp.exp(p*mp.log(abs(mp.cos(th - i*mp.pi/m)))) for i in range(m))
def r_and_rp(th, p):
    f = fsum(th,p)
    fp = sum(-mp.sin(th - i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    r = f**(-1/p)
    rp = -(1/p)*f**(-1/p-1)*fp
    return r, rp
# Integrand over the chamber centered at 0: theta in [-pi/14, pi/14]
def integrand(th, p):
    r, rp = r_and_rp(th,p)
    return mp.sqrt(r*r+rp*rp)
# Perimeter via tanh-sinh with explicit high precision; full chamber width pi/7
def perim(p, deg=8):
    a = -mp.pi/14; b = mp.pi/14
    val = mp.quad(lambda th: integrand(th,p), [a,b], method='tanh-sinh', maxdegree=deg)
    return 2*m*val
P1 = perim(1)
print("P(1) via quad =", mp.nstr(P1,40))
print("P(1)-14L      =", mp.nstr(P1-14*LL,10))
open('progress.log','a').write('PG6 s_A1: P(1) quad vs 14L diff = %s\n'%mp.nstr(abs(P1-14*LL),3))
