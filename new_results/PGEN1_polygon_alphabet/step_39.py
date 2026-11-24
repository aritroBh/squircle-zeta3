# step_16: regular 2m-gon edge geometry: kappa_j, marked points p_j = -1/kappa_j, for m=4,5,6
import sympy as sp

def regular_polygon_edges(m):
    """K_p = {sum_{j=1}^m |<n_j,x>|^p <= 1}, n_j = (cos(j pi/m), sin(j pi/m)).
    Return corners (exact algebraic) and per-edge term data for one edge."""
    pi = sp.pi
    n = [sp.Matrix([sp.cos(j*pi/m), sp.sin(j*pi/m)]) for j in range(1, m+1)]
    # corner on kink line j: x = R*(-sin th_j, cos th_j), R = 1/sum_{i != j} |<n_i, dir>|
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
    # edge from corner 0 (kink line n_1) to corner 1 (kink line n_2)
    V0, V1 = corners[0], corners[1]
    print(f"--- regular {2m}-gon: edge from kink line 1 to kink line 2")
    data = []
    for i in range(m):
        a0 = sp.simplify(abs(n[i].dot(V0)))
        a1 = sp.simplify(abs(n[i].dot(V1)))
        # sign of <n_i, y(t)> constant on edge
        if a0 == 0:
            data.append((i+1, 'corner0', None, None))
        elif a1 == 0:
            data.append((i+1, 'corner1', None, None))
        else:
            r = sp.simplify(a1/a0)
            kap = sp.simplify(r-1)
            data.append((i+1, 'far', kap, r))
    for d in data:
        print("   n_%d %s  kappa=%s r=%s" % d)
