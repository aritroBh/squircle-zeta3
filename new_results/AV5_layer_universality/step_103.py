import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 70
def spow(t, a):
    if t == 0: return mpf(0)
    return (1 if t>0 else -1)*abs(t)**a
def h_oct(c1, c2, p):
    q = p/(p-1)
    sq2 = sqrt(2)
    def eqs(z3, z4):
        z1 = c1 - (z3+z4)/sq2
        z2 = c2 - (z3-z4)/sq2
        f1 = -(spow(z1,q-1)+spow(z2,q-1))/sq2 + spow(z3,q-1)
        f2 = -(spow(z1,q-1)-spow(z2,q-1))/sq2 + spow(z4,q-1)
        return (f1, f2)
    z3, z4 = c1/3, c2/3
    for _ in range(60):
        f1, f2 = eqs(z3, z4)
        # Newton with numeric Jacobian
        h = mpf('1e-30')
        f1a, f2a = eqs(z3+h, z4); f1b, f2b = eqs(z3, z4+h)
        J11, J12 = (f1a-f1)/h, (f1b-f1)/h
        J21, J22 = (f2a-f1)/h, (f2b-f2)/h
        det = J11*J22 - J12*J21
        dz3 = (-f1*J22 + f2*J12)/det
        dz4 = (-J11*f2 + J21*f1)/det
        z3 += dz3; z4 += dz4
        if abs(dz3)+abs(dz4) < mpf(10)**(-55): break
    z1 = c1 - (z3+z4)/sq2; z2 = c2 - (z3-z4)/sq2
    return (abs(z1)**q+abs(z2)**q+abs(z3)**q+abs(z4)**q)**(1/q)
# validate p=2: circle radius 1/sqrt2: h = 1/sqrt2, P = 2 pi /sqrt2
for th in [mpf(0), mpf('0.6'), mpf('1.2')]:
    print(nstr(h_oct(cos(th),sin(th),mpf(2)),25), " should be 1/sqrt2 =", nstr(1/sqrt(2),25))
Pfull = mpmath.quad(lambda th: h_oct(cos(th),sin(th),mpf(2)), [0, pi/2, pi, 3*pi/2, 2*pi])
print("P(2) =", nstr(Pfull,30), " pi*sqrt2:", nstr(pi*sqrt(2),30))
# octagon limit perimeter: vertices (±1,±(√2-1)),(±(√2-1),±1): 8 edges of length 2(√2-1)
print("oct limit:", nstr(8*2*(sqrt(2)-1),20))
