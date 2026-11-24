# Independent referee check: W4 from its definition vs R transcribed by hand from main.tex (T3).
import mpmath as mp
mp.mp.dps = 90
r = mp.sqrt(2); a = 1+r; b = r-1; L = mp.log(a); lam = mp.log(2); pi = mp.pi
Li = mp.polylog; z3 = mp.zeta(3)
def l1(v): return mp.log((v-b)*(v+a)/(4*v))
def l2(v): return mp.log((a-v)*(v+b)/(4*v))
def f0(v): return l1(v)**3*(mp.mpf(1)/(6*v) - 1/(3*(v+a)))          # rational prefactor part
def f1(v): return -mp.mpf(7)/48*l1(v)**3/v + l1(v)**2*l2(v)/(16*v)  # sqrt2 prefactor part
pts = [b, (b+1)/2, 1, (1+a)/2, a]
U = mp.quad(f0, pts, method='tanh-sinh'); V = mp.quad(f1, pts, method='tanh-sinh')
W4 = U + r*V
def Rparts():
    # returns list of (p, q, value): coefficient p + q*sqrt2
    F = mp.mpf
    T = []
    T += [(1,-573,Li(4,b)), (-1,-573,Li(4,-b)), (0,56,Li(4,-b**2)), (0,17,Li(4,1/r)),
          (F(-37)/4,F(-17)/16,Li(4,F(1)/2)), (10,86,Li(4,1-1/r)), (-10,86,Li(4,2-r))]
    T += [(40,7,L*Li(3,b)), (-40,7,L*Li(3,-b)), (-72,0,L*Li(3,1/r)), (F(63)/8,F(-3)/16,L*z3)]
    T += [(6,13,lam*Li(3,b)), (6,-13,lam*Li(3,-b)), (0,F(-47)/2,lam*Li(3,1/r)), (F(-3)/8,F(329)/128,lam*z3)]
    T += [(30,0,Li(2,b)**2)]
    T += [(30,F(7)/2,Li(2,b)*L**2), (-24,0,Li(2,b)*lam*L), (0,F(-19)/8,Li(2,b)*lam**2), (F(-19)/4,F(1)/24,Li(2,b)*pi**2)]
    T += [(F(103)/384,F(-4963)/5760,pi**4)]
    T += [(F(67)/24,F(185)/96,pi**2*L**2), (F(-35)/12,F(41)/12,pi**2*lam*L), (F(5)/12,F(-1105)/384,pi**2*lam**2)]
    T += [(F(13)/24,F(-25)/24,L**4), (F(-73)/6,F(-29)/3,L**3*lam), (F(63)/8,F(67)/16,L**2*lam**2),
          (F(41)/24,F(17)/6,L*lam**3), (F(-61)/192,F(137)/192,lam**4)]
    return T
T = Rparts()
R0 = sum(mp.mpf(p)*v for p,q,v in T); R1 = sum(mp.mpf(q)*v for p,q,v in T)
R = R0 + r*R1
print('W4 =', mp.nstr(W4, 70)); print('R  =', mp.nstr(R, 70))
print('W4-R =', mp.nstr(W4-R, 5), ' U-R0 =', mp.nstr(U-R0,5), ' V-R1 =', mp.nstr(V-R1,5))
print('U =', mp.nstr(U,50)); print('V =', mp.nstr(V,50))
for d in (W4-R, U-R0, V-R1):
    print('digits agreement ~', mp.nstr(-mp.log10(abs(d)) if d else mp.inf, 4))
