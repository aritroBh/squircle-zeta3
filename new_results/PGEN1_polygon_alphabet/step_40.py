# step_16 (fixed): regular 2m-gon edge geometry
import sympy as sp

def regular_polygon_edges(m):
    pi = sp.pi
    n = [sp.Matrix([sp.cos(j*pi/m), sp.sin(j*pi/m)]) for j in range(1, m+1)]
    corners = []
    for j in range(m):
        th = (j+1)*pi/m
        d = sp.Matrix([-sp.sin(th), sp.cos(th)])
        s = 0
        for i in range(m):
            s += sp.Abs(n[i].dot(d))
        R = 1/s
        corners.append(sp.simplify(R*d))
    return n, corners

for m in (4,5,6):
    n, corners = regular_polygon_edges(m)
    V0, V1 = corners[0], corners[1]
    print("--- regular %d-gon: edge from kink line 1 to kink line 2" % (2*m))
    for i in range(m):
        a0 = sp.simplify(abs(n[i].dot(V0)))
        a1 = sp.simplify(abs(n[i].dot(V1)))
        if a0 == 0:
            print("   n_%d corner-at-V0" % (i+1))
        elif a1 == 0:
            print("   n_%d corner-at-V1" % (i+1))
        else:
            r = sp.simplify(a1/a0)
            kap = sp.simplify(r-1)
            p = sp.simplify(-1/kap)
            print("   n_%d far: kappa=%s r=%s p=-1/kappa=%s" % (i+1, kap, r, p))
