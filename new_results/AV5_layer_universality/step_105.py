import mpmath, json, time
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr
mp.dps = 70
def spow(t, a):
    if t == 0: return mpf(0)
    return (1 if t>0 else -1)*abs(t)**a
def h_oct(c1, c2, p):
    q = p/(p-1); sq2 = sqrt(2)
    def G(z3, z4):
        z1 = c1 - (z3+z4)/sq2
        z2 = c2 - (z3-z4)/sq2
        return abs(z1)**q+abs(z2)**q+abs(z3)**q+abs(z4)**q
    def nested():
        # minimize G: nested bisection
        z4lo, z4hi = -3*max(abs(c1),abs(c2)), 3*max(abs(c1),abs(c2))
        def min_z3(z4):
            lo, hi = -3*max(abs(c1),abs(c2)), 3*max(abs(c1),abs(c2))
            gr = (sqrt(5)-1)/2
            x1, x2 = hi-gr*(hi-lo), lo+gr*(hi-lo)
            f1, f2 = G(x1,z4), G(x2,z4)
            for _ in range(120):
                if f1 < f2:
                    lo, x1, f1 = x1, x2, f2
                    x2 = lo+gr*(hi-lo); f2 = G(x2,z4)
                else:
                    hi, x2, f2 = x2, x1, f1
                    x1 = hi-gr*(hi-lo); f1 = G(x1,z4)
            return (f1+f2)/2
        gr = (sqrt(5)-1)/2
        lo, hi = z4lo, z4hi
        x1, x2 = hi-gr*(hi-lo), lo+gr*(hi-lo)
        f1, f2 = min_z3(x1), min_z3(x2)
        for _ in range(60):
            if f1 < f2:
                lo, x1, f1 = x1, x2, f2
                x2 = lo+gr*(hi-lo); f2 = min_z3(x2)
            else:
                hi, x2, f2 = x2, x1, f1
                x1 = hi-gr*(hi-lo); f1 = min_z3(x1)
        return (f1+f2)/2
    # try Newton from several guesses, else nested
    for g3, g4 in [(c1/3,c2/3),(c1/2,c2/2),(mpf(0),mpf(0)),(c1,c2)]:
        try:
            z3, z4 = mpf(g3), mpf(g4)
            hh = mpf('1e-25')
            ok = False
            for _ in range(50):
                z1 = c1 - (z3+z4)/sq2; z2 = c2 - (z3-z4)/sq2
                f1 = -(spow(z1,q-1)+spow(z2,q-1))/sq2 + spow(z3,q-1)
                f2 = -(spow(z1,q-1)-spow(z2,q-1))/sq2 + spow(z4,q-1)
                if abs(f1)+abs(f2) < mpf(10)**(-50): ok=True; break
                def eqs(a,b):
                    zz1 = c1 - (a+b)/sq2; zz2 = c2 - (a-b)/sq2
                    return (-(spow(zz1,q-1)+spow(zz2,q-1))/sq2 + spow(a,q-1),
                            -(spow(zz1,q-1)-spow(zz2,q-1))/sq2 + spow(b,q-1))
                f1a, f2a = eqs(z3+hh, z4); f1b, f2b = eqs(z3, z4+hh)
                J11, J12 = (f1a-f1)/hh, (f1b-f1)/hh
                J21, J22 = (f2a-f1)/hh, (f2b-f2)/hh
                det = J11*J22 - J12*J21
                if det == 0: break
                dz3 = (-f1*J22 + f2*J12)/det; dz4 = (-J11*f2 + J21*f1)/det
                # damping
                mx = max(abs(dz3),abs(dz4),1)
                if mx > 1: dz3, dz4 = dz3/mx, dz4/mx
                z3 += dz3; z4 += dz4
            if ok:
                return G(z3,z4)**(1/q)
        except Exception:
            pass
    return nested()**(1/q)
# validate
print("p=2 circle check:", nstr(h_oct(cos(mpf('0.6')),sin(mpf('0.6')),mpf(2)),25), nstr(1/sqrt(2),25))
Pfull = mpmath.quad(lambda th: h_oct(cos(th),sin(th),mpf(2)), [0, pi/2, pi, 3*pi/2, 2*pi])
print("P(2):", nstr(Pfull,25), nstr(pi*sqrt(2),25))
