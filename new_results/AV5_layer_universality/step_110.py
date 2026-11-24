import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 50
def spow(t, a):
    if t == 0: return mpf(0)
    return (1 if t>0 else -1)*abs(t)**a
def h_oct(c1, c2, p):
    q = p/(p-1); sq2 = sqrt(2)
    def spowd(t, a):
        if t == 0: return mpf(0)
        return a*abs(t)**(a-1)
    def solve(g3, g4):
        z3, z4 = mpf(g3), mpf(g4)
        for _ in range(80):
            z1 = c1 - (z3+z4)/sq2; z2 = c2 - (z3-z4)/sq2
            f1 = -(spow(z1,q-1)+spow(z2,q-1))/sq2 + spow(z3,q-1)
            f2 = -(spow(z1,q-1)-spow(z2,q-1))/sq2 + spow(z4,q-1)
            if abs(f1)+abs(f2) < mpf(10)**(-40): return (z3, z4)
            a1, a2 = spowd(z1,q-1), spowd(z2,q-1)
            J11 = (a1+a2)/2 + spowd(z3,q-1)
            J22 = (a1+a2)/2 + spowd(z4,q-1)
            J12 = (a1-a2)/2; J21 = (a1-a2)/2
            det = J11*J22 - J12*J21
            if det == 0: return None
            dz3 = (-f1*J22 + f2*J12)/det; dz4 = (-J11*f2 + J21*f1)/det
            mx = max(abs(dz3), abs(dz4))
            if mx > mpf('0.5'): dz3, dz4 = dz3*mpf('0.5')/mx, dz4*mpf('0.5')/mx
            z3 += dz3; z4 += dz4
        return None
    r = None
    for g in [(c1/3,c2/3),(c1/2,c2/2),(0,0),(c1,c2),(-c1,-c2),(c1/4,-c2/4)]:
        r = solve(*g)
        if r: break
    if not r: raise RuntimeError("no converge")
    z3, z4 = r
    z1 = c1 - (z3+z4)/sq2; z2 = c2 - (z3-z4)/sq2
    return (abs(z1)**q+abs(z2)**q+abs(z3)**q+abs(z4)**q)**(1/q)
print("p=2:", nstr(h_oct(cos(mpf('0.6')),sin(mpf('0.6')),mpf(2)),25), nstr(1/sqrt(2),25))
def P_oct(p):
    return 16*mpmath.quad(lambda th: h_oct(cos(th),sin(th),p), [0, pi/16, pi/8])
t0=time.time()
print("P(2):", nstr(P_oct(mpf(2)),35), nstr(pi*sqrt(2),35), "time", time.time()-t0)
