# s04: exact integration of the LOCAL part of int_0^1 E_k dt (all non-vanishing logs frozen at their
# constant values), via the Beta table.  For polygons in which every edge is parallel to all the other
# kink lines it does not touch (rhombus; affine-regular hexagon), the local part IS the whole edge
# integrand, so this gives P^{(k)}(1) in closed form.
import sympy as sp, pickle, sys
from s03_beta_table import J, z2, z3

D = pickle.load(open('s01_E.pkl', 'rb'))
E = D['E']
t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
l0, l1 = sp.symbols('l0 l1')
tau = sp.symbols('tau1:4'); gam = sp.symbols('gam1:4')
X, Y = sp.symbols('X Y')  # X = log t, Y = log(1-t)

def edge_exact(V0, V1, ns, ws, lams):
    """exact edge data; returns dict and list of non-local normals (index, a_i(t), a_i', mu_i, r_i)"""
    c0 = [sp.nsimplify(sp.simplify(n[0]*V0[0] + n[1]*V0[1])) for n in ns]
    c1 = [sp.nsimplify(sp.simplify(n[0]*V1[0] + n[1]*V1[1])) for n in ns]
    j0 = [i for i in range(len(ns)) if sp.simplify(c0[i]) == 0]; j1 = [i for i in range(len(ns)) if sp.simplify(c1[i]) == 0]
    assert len(j0) == 1 and len(j1) == 1, (c0, c1)
    j0, j1 = j0[0], j1[0]
    Dv = (V1[0] - V0[0], V1[1] - V0[1])
    Le = sp.sqrt(sp.simplify(Dv[0]**2 + Dv[1]**2))
    tv = (Dv[0]/Le, Dv[1]/Le)
    he = sp.simplify(sp.Abs(V0[0]*tv[1] - V0[1]*tv[0]))
    s0e = sp.simplify(V0[0]*tv[0] + V0[1]*tv[1])
    dat = dict(alpha=sp.simplify(ws[j0]*sp.Abs(c1[j0])), beta=sp.simplify(ws[j1]*sp.Abs(c0[j1])), L=Le, h=he, s0=s0e,
               mu0=sp.log(sp.Abs(c1[j0])) + lams[j0], mu1=sp.log(sp.Abs(c0[j1])) + lams[j1])
    nl = []
    for i in range(len(ns)):
        if i in (j0, j1): continue
        sg = sp.sign(c0[i])
        a_t = sp.simplify(ws[i]*sg*(c0[i] + t*(c1[i] - c0[i])))
        ap = sp.simplify(ws[i]*sg*(c1[i] - c0[i]))
        nl.append(dict(i=i, a=a_t, ap=ap, mu=sp.log(sp.Abs(c0[i])) + lams[i], r=sp.simplify(c1[i]/c0[i])))
    return dat, nl

def local_integral(dat, nl, k):
    """exact int_0^1 of E_k with non-local logs l_i -> mu_i (their value; exact if r_i = 1)"""
    sub = {al: dat['alpha'], be: dat['beta'], L: dat['L'], h: dat['h'], s0: dat['s0'],
           l0: X + dat['mu0'], l1: Y + dat['mu1']}
    for kk in range(1, 4):
        sub[tau[kk-1]] = sum(q['a']*q['mu']**kk for q in nl)
        sub[gam[kk-1]] = sum(q['ap']*q['mu']**kk for q in nl)
    ex = sp.expand(E[k].xreplace(sub))
    P = sp.Poly(ex, t, X, Y)
    tot = 0
    for (n, a, b), c in P.terms():
        tot += c*J(n, 0, a, b)
    return sp.expand(tot)

def polygon_exact(ns, ws, lams, dirs):
    """dirs: exact unit vertex directions in ccw order covering the full turn"""
    V = []
    for dv in dirs:
        S = sum(w*sp.Abs(n[0]*dv[0] + n[1]*dv[1]) for n, w in zip(ns, ws))
        rho = sp.nsimplify(sp.simplify(1/S))
        V.append((sp.simplify(rho*dv[0]), sp.simplify(rho*dv[1])))
    return V

def closed_form(ns, ws, lams, dirs, name, kmax=3):
    V = polygon_exact(ns, ws, lams, dirs)
    res = [0]*(kmax+1); z3part = 0; allr1 = True; Zlaw = 0
    for q in range(len(V)):
        dat, nl = edge_exact(V[q], V[(q+1) % len(V)], ns, ws, lams)
        for nn in nl:
            if sp.simplify(nn['r'] - 1) != 0: allr1 = False
        for k in range(kmax+1):
            res[k] += local_integral(dat, nl, k)
        A, B, Le, he = dat['alpha'], dat['beta'], dat['L'], dat['h']
        u0 = -dat['s0']/Le; u1 = (dat['s0'] + Le)/Le
        Zlaw += 6*A*B*he**2/Le*(3*(A*u0 + B*u1) - 2)
    out = [sp.simplify(sp.factorial(k)*res[k]) for k in range(kmax+1)]
    print(f"== {name}: vertices {V}")
    print("   all non-local ratios r_i = 1 (fully local, closed form exact):", allr1)
    for k in range(kmax+1):
        print(f"   P^({k})(1) =", sp.collect(sp.expand(out[k]), [z3, z2]))
    print("   z3-coefficient of P'''(1):", sp.nsimplify(sp.simplify(sp.expand(out[3]).coeff(z3))),
          "  edge-law prediction:", sp.nsimplify(sp.simplify(Zlaw)),
          "  difference:", sp.simplify(sp.expand(out[3]).coeff(z3) - Zlaw))
    return out, allr1, sp.simplify(Zlaw)

if __name__ == '__main__':
    r2, r3 = sp.sqrt(2), sp.sqrt(3)
    results = {}
    # rhombus |x/a|^p + |y/b|^p <= 1 (power convention): n=(0,1) weight 1/b lam=-log b; n=(1,0) weight 1/a lam=-log a
    a_, b_ = sp.Integer(1), sp.Integer(2)
    ns = [(0, 1), (1, 0)]; ws = [1/b_, 1/a_]; lams = [-sp.log(b_), -sp.log(a_)]
    dirs = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    out, ok, Z = closed_form(ns, ws, lams, dirs, "rhombus a=1,b=2 power convention")
    r = sp.sqrt(a_**2 + b_**2)
    print("   R2 target P'''(1) = 24a^2b^2/r^3 (z3 - pi^2/3):  difference =",
          sp.simplify(out[3].subs(z2, sp.pi**2/6) - 24*a_**2*b_**2/r**3*(z3 - sp.pi**2/3)))
    results['rhombus12_power'] = [str(x) for x in out]
    # same rhombus, problem convention (weights fixed): |y|^p/b + |x|^p/a <= 1
    out, ok, Z = closed_form(ns, ws, [0, 0], dirs, "rhombus a=1,b=2 fixed-weight convention")
    results['rhombus12_fixed'] = [str(x) for x in out]
    # regular hexagon, m = 3, problem convention: n_i = (-sin b_i, cos b_i), b_i = 0, -pi/3, -2pi/3
    ns = [(0, 1), (r3/2, sp.Rational(1, 2)), (r3/2, -sp.Rational(1, 2))]
    ws = [1, 1, 1]; lams = [0, 0, 0]
    dirs = [(sp.cos(k*sp.pi/3), sp.sin(k*sp.pi/3)) for k in range(6)]
    out, ok, Z = closed_form(ns, ws, lams, dirs, "regular hexagon m=3")
    l2 = sp.log(2)
    D15 = r3/8*(16*l2**3 + 3*sp.pi**2*l2 - 96*l2**2 + 96*l2 - 9*z3 - 6*sp.pi**2)
    print("   D15 PSLQ formula difference:", sp.simplify(sp.expand(out[3].subs(z2, sp.pi**2/6) - D15)))
    results['hex_regular'] = [str(x) for x in out]
    # affine-regular hexagon H_1 = {|x|+|y|+|x-y| <= 1}; kinks y=0 (n=(0,1)), x=0 (n=(1,0)), x=y (n=(-1,1)/sqrt2)
    ns = [(0, 1), (1, 0), (-1/r2, 1/r2)]
    dirs = [(1, 0), (1/r2, 1/r2), (0, 1), (-1, 0), (-1/r2, -1/r2), (0, -1)]
    out, ok, Z = closed_form(ns, [1, 1, r2], [0, 0, 0], dirs, "H_p fixed-weight: |x|^p+|y|^p+sqrt2|(x-y)/sqrt2|^p<=1")
    results['H_fixed'] = [str(x) for x in out]
    out, ok, Z = closed_form(ns, [1, 1, r2], [0, 0, sp.log(r2)], dirs, "H_p power: |x|^p+|y|^p+|x-y|^p<=1")
    results['H_power'] = [str(x) for x in out]
    pickle.dump(results, open('s04_closed_forms.pkl', 'wb'))
    with open('s04_closed_forms.txt', 'w') as f:
        for k, v in results.items():
            f.write(k + '\n' + '\n'.join(f'  P^({i})(1) = {x}' for i, x in enumerate(v)) + '\n')
    print("saved s04_closed_forms.pkl/.txt")
