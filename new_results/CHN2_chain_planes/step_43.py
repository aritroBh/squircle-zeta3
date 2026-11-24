
import mpmath as mp
mp.mp.dps = 70
EPS = mp.mpf(10)**(-60)
def ts1d(f, h=mp.mpf('0.08')):
    s = mp.mpf(0)
    for k in range(0, 600):
        t = k*h; u_ = mp.mpf(mp.pi)/2*mp.sinh(t); y = mp.tanh(u_)
        if y >= 1: break
        w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
        term = f(y)*w; s += term
        if k>25 and abs(term) < EPS: break
    for k in range(1, 600):
        t = -k*h; u_ = mp.mpf(mp.pi)/2*mp.sinh(t); y = mp.tanh(u_)
        if y <= -1: break
        w = (mp.mpf(mp.pi)/2*mp.cosh(t))/(mp.cosh(u_)**2)
        term = f(y)*w; s += term
        if k>25 and abs(term) < EPS: break
    return s*h
def int_ab(f,a,b): return ts1d(lambda y: f(a+(y+1)/2*(b-a)))*(b-a)/2
def I2s(T,a1,a2):
    return int_ab(lambda w: (1-w**(1/a2))**a1, 0, T**a2)
def vol_t(u1,u2,u3,u4,u5):
    f = lambda s: I2s(1-s**(1/u3),u1,u2)*I2s(1-s**(1/u3),u5,u4)
    return int_ab(f, 0, mp.mpf('0.5')) + int_ab(f, mp.mpf('0.5'), 1)

# direction (0,1,1,1,1)*u ; I2 with a1=0: int_0^T a2 t^(a2-1) dt = T^a2 -- fine
a=(mp.mpf(0),mp.mpf(1),mp.mpf(1),mp.mpf(1),mp.mpf(1))
hs=[mp.mpf(1)/100, mp.mpf(1)/200, mp.mpf(1)/400, mp.mpf(1)/800, mp.mpf(1)/1600]
Vs=[vol_t(*[x*h for x in a]) for h in hs]
# V = 1 + c2 u^2 + c3 u^3 + c4 u^4 + c5 u^5 + c6 u^6
Amat = mp.matrix(5,5)
for i,h in enumerate(hs):
    for j in range(5):
        Amat[i,j] = h**(j+2)
b = mp.matrix([V-1 for V in Vs])
sol = mp.lu_solve(Amat, b)
z3 = mp.zeta(3)
for j,c in enumerate(sol):
    print(f"c{j+2} = {mp.nstr(c,25)}   /zeta({j+2}) = {mp.nstr(c/mp.zeta(j+2),15)}" if j+2 in (2,3,4) else f"c{j+2}={mp.nstr(c,20)}")
print("c3/zeta3:", mp.nstr(sol[1]/z3, 20), " -> mine said 12, CHN1 said 10")
