import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 70
def spow(t, a):
    if t == 0: return mpf(0)
    return (1 if t>0 else -1)*abs(t)**a
def h_oct(c1, c2, p):
    q = p/(p-1); sq2 = sqrt(2)
    def eqs(z3, z4):
        z1 = c1 - (z3+z4)/sq2
        z2 = c2 - (z3-z4)/sq2
        return (-(spow(z1,q-1)+spow(z2,q-1))/sq2 + spow(z3,q-1),
                -(spow(z1,q-1)-spow(z2,q-1))/sq2 + spow(z4,q-1))
    z3, z4 = c1/3, c2/3
    h = mpf('1e-30')
    for _ in range(60):
        f1, f2 = eqs(z3, z4)
        f1a, f2a = eqs(z3+h, z4); f1b, f2b = eqs(z3, z4+h)
        J11, J12 = (f1a-f1)/h, (f1b-f1)/h
        J21, J22 = (f2a-f1)/h, (f2b-f2)/h
        det = J11*J22 - J12*J21
        dz3 = (-f1*J22 + f2*J12)/det; dz4 = (-J11*f2 + J21*f1)/det
        z3 += dz3; z4 += dz4
        if abs(dz3)+abs(dz4) < mpf(10)**(-55): break
    z1 = c1 - (z3+z4)/sq2; z2 = c2 - (z3-z4)/sq2
    return (abs(z1)**q+abs(z2)**q+abs(z3)**q+abs(z4)**q)**(1/q)
def P_oct(p):
    return mpmath.quad(lambda th: h_oct(cos(th),sin(th),p),
                       [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 2*pi])
t0=time.time()
epss = [mpf('0.02')/mpf('1.6')**j for j in range(7)]
Ps = [P_oct(1/e) for e in epss]
json.dump([nstr(e,30) for e in epss], open('ckpt_oct_epss.json','w'))
json.dump([nstr(P,65) for P in Ps], open('ckpt_oct_Ps.json','w'))
Poct = 16*(sqrt(2)-1)
for e,P in zip(epss,Ps):
    print(f"eps={nstr(e,10)} Poct-P={nstr(Poct-P,12)}")
print("time", time.time()-t0)
