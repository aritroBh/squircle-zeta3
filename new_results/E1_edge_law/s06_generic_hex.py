# s06: generic (non-affine-regular) hexagons.  Decompose 6*sum_e int_0^1 E_3 dt into
#   (i)  local part (non-local log frozen): exact via Beta table = Q0 + Q2 z2 + Q3 z3, with Q3 = edge law;
#   (ii) non-local part (terms with >= 1 power of N(t) = log(1+(r-1)t)): by quadrature AND as a power series
#        in x = r-1 (or x' = 1/r-1) whose coefficients use only Beta moments with a+b <= 2 (no zeta(3) possible).
# Compare the total with D15b's polar coefficient integrals (independent derivation) at high precision.
import sys, json, pickle, mpmath as mp, sympy as sp
sys.path.insert(0, 'research/08_theory/D15b_polygon_law')
import poly_coeffs as pc
from edge_num import edges, zeta3_edge_law
from s03_beta_table import J, z2, z3

mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 50
D = pickle.load(open('s01_E.pkl', 'rb')); E3 = D['E'][3]
t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
l0, l1 = sp.symbols('l0 l1'); tau = sp.symbols('tau1:4'); gam = sp.symbols('gam1:4')
X, Y, Nn, mu0, mu1, mu2 = sp.symbols('X Y N mu0 mu1 mu2')
# m = 3: one non-local normal per edge, a_2(t) = 1 - beta + (beta - alpha) t, a_2' = beta - alpha
a2 = 1 - be + (be - al)*t; a2p = be - al
sub = {l0: X + mu0, l1: Y + mu1}
for k in range(1, 4):
    sub[tau[k-1]] = a2*(mu2 + Nn)**k; sub[gam[k-1]] = a2p*(mu2 + Nn)**k
poly = sp.Poly(sp.expand(E3.xreplace(sub)), t, X, Y, Nn)
terms = [(mon, sp.lambdify((al, be, L, h, s0, mu0, mu1, mu2), c, 'mpmath')) for mon, c in poly.terms()]
print("E3 for m=3 edge: number of monomials t^n X^a Y^b N^j:", len(terms),
      " max j:", max(m[3] for m, _ in terms), " max a+b+j:", max(m[1]+m[2]+m[3] for m, _ in terms))
assert all(m[1] + m[2] + m[3] <= 3 for m, _ in terms)

HC = {1: [mp.mpf(0)], 2: [mp.mpf(0)]}
def Hn(q, r):
    arr = HC[r]
    while len(arr) <= q:
        arr.append(arr[-1] + mp.mpf(1)/len(arr)**r)
    return arr[q]
def Jnum2(n, k, a, b):
    """int t^n (1-t)^k log^a t log^b(1-t), a+b <= 2: only psi, psi' needed -> rationals and zeta(2)"""
    assert a + b <= 2
    B = mp.factorial(n)*mp.factorial(k)/mp.factorial(n+k+1)
    fx = Hn(n, 1) - Hn(n+k+1, 1); fy = Hn(k, 1) - Hn(n+k+1, 1)
    fxx = -Hn(n, 2) + Hn(n+k+1, 2); fyy = -Hn(k, 2) + Hn(n+k+1, 2)
    fxy = -(mp.zeta(2) - Hn(n+k+1, 2))
    return {(0, 0): B, (1, 0): B*fx, (0, 1): B*fy, (2, 0): B*(fx**2 + fxx),
            (0, 2): B*(fy**2 + fyy), (1, 1): B*(fx*fy + fxy)}[(a, b)]

def series_nonlocal(coef, r):
    """sum over monomials with j>=1 of c * int t^n X^a Y^b N^j, N = log(1+(r-1)t), by power series"""
    x = r - 1
    rho = min(abs(x), abs(1/r - 1))
    nmax = int((mp.mp.dps + 10)/(-mp.log10(rho))) + 30
    if abs(x) <= abs(1/r - 1):
        # N = sum_{q>=1} (-1)^{q+1} x^q t^q / q ; N^j as series in t with coefficients in x
        base = [mp.mpf(0)] + [(-1)**(q+1)*x**q/q for q in range(1, nmax+1)]
        mode = 't'
    else:
        # N = log r + log(1 + (1/r - 1)(1-t)) -> series in s = 1-t
        xp = 1/r - 1
        base = [mp.log(r)] + [(-1)**(q+1)*xp**q/q for q in range(1, nmax+1)]
        mode = 's'
    pw = {0: [mp.mpf(1)] + [mp.mpf(0)]*nmax}
    for j in range(1, 4):
        prev = pw[j-1]
        pw[j] = [mp.fsum([prev[i]*base[q-i] for i in range(q+1)]) for q in range(nmax+1)]
    tot = mp.mpf(0)
    for (n, a, b, j), c in coef:
        if j == 0: continue
        for q in range(nmax+1):
            cq = pw[j][q]
            if cq == 0: continue
            tot += c*cq*(Jnum2(n + q, 0, a, b) if mode == 't' else Jnum2(n, q, a, b))
    return tot, mode

def run(name, bs, ws, lams):
    Z, parts = zeta3_edge_law(bs, ws)
    Es = edges(bs, ws, lams)
    loc = 0; nl_quad = mp.mpf(0); nl_ser = mp.mpf(0); full_quad = mp.mpf(0); rs = []
    for e in Es:
        j0, j1 = e['j0'], e['j1']; i2 = [i for i in range(3) if i not in (j0, j1)][0]
        A, B = e['alpha'], e['beta']
        m0 = mp.log(abs(e['c1'][j0])) + lams[j0]; m1 = mp.log(abs(e['c0'][j1])) + lams[j1]
        m2 = mp.log(abs(e['c0'][i2])) + lams[i2]; r = e['c1'][i2]/e['c0'][i2]
        assert abs(ws[i2]*abs(e['c0'][i2]) - (1 - B)) < mp.mpf(10)**(-mp.mp.dps+5)
        assert abs(r - (1 - A)/(1 - B)) < mp.mpf(10)**(-mp.mp.dps+5)
        rs.append(r)
        coef = [(mon, f(A, B, e['L'], e['h'], e['s0'], m0, m1, m2)) for mon, f in terms]
        # local exact (sympy exact Beta values, evaluated numerically at the end)
        lq0 = mp.mpf(0); lq2 = mp.mpf(0); lq3 = mp.mpf(0)
        for (n, a, b, j), c in coef:
            if j: continue
            ex = J(n, 0, a, b)
            lq0 += c*mp.mpf(str(sp.N(ex.subs({z2: 0, z3: 0}), mp.mp.dps + 10)))
            lq2 += c*mp.mpf(str(sp.N(ex.coeff(z2), mp.mp.dps + 10)))
            lq3 += c*mp.mpf(str(sp.N(ex.coeff(z3), mp.mp.dps + 10)))
        Nf = lambda tt: mp.log(1 + (r - 1)*tt)
        def integrand(tt, only_nl):
            Xv, Yv, Nv = mp.log(tt), mp.log(1 - tt), Nf(tt)
            return mp.fsum([c*tt**n*Xv**a*Yv**b*Nv**j for (n, a, b, j), c in coef if (j > 0 or not only_nl)])
        nlq = mp.quad(lambda tt: integrand(tt, True), [0, mp.mpf(1)/2, 1])
        fq = mp.quad(lambda tt: integrand(tt, False), [0, mp.mpf(1)/2, 1])
        nls, mode = series_nonlocal(coef, r)
        e.update(loc=(lq0, lq2, lq3), nlq=nlq, nls=nls, mode=mode, fq=fq)
        nl_quad += nlq; nl_ser += nls; full_quad += fq
        print(f"  edge {j0}->{j1}: alpha={mp.nstr(A,8)} beta={mp.nstr(B,8)} r={mp.nstr(r,8)} [{mode}]  "
              f"z3-coef(local)*6 = {mp.nstr(6*lq3, 20)}  edge law Z_e = {mp.nstr(parts[Es.index(e)], 20)}  "
              f"|nl quad - nl series| = {mp.nstr(abs(nlq - nls), 3)}")
    Q0 = 6*mp.fsum([e['loc'][0] for e in Es]); Q2 = 6*mp.fsum([e['loc'][1] for e in Es])
    Q3 = 6*mp.fsum([e['loc'][2] for e in Es])
    P3_dec = Q0 + Q2*mp.zeta(2) + Q3*mp.zeta(3) + 6*nl_ser
    # independent: D15b polar coefficient integrals
    pcv = pc.P_coeffs_general(bs, ws, lams, kmax=3)
    P3_polar = 6*pcv[3]
    print(f" {name}: Z (edge law) = {mp.nstr(Z, 30)};  6*sum local z3-coef = {mp.nstr(Q3, 30)};  diff {mp.nstr(Z - Q3, 3)}")
    print(f"   P'''(1) decomposed  = {mp.nstr(P3_dec, 40)}")
    print(f"   P'''(1) polar D15b  = {mp.nstr(P3_polar, 40)}   |diff| = {mp.nstr(abs(P3_dec - P3_polar), 3)}")
    print(f"   remainder R = P''' - Z zeta(3) = Q0 + Q2 zeta(2) + 6*NL, with Q0 = {mp.nstr(Q0, 25)}, Q2 = {mp.nstr(Q2, 25)}, 6*NL = {mp.nstr(6*nl_ser, 25)}")
    return dict(Z=str(Z), Q0=str(Q0), Q2=str(Q2), Q3=str(Q3), NL=str(6*nl_ser), P3_dec=str(P3_dec),
                P3_polar=str(P3_polar), r=[str(x) for x in rs], parts=[str(x) for x in parts])

if __name__ == '__main__':
    out = {}
    out['hexA_w123'] = run('hexA_w123', [mp.mpf(0), -mp.pi/3, -2*mp.pi/3], [1, 2, 3], [0, 0, 0])
    out['hexB_angles'] = run('hexB_angles', [mp.mpf(0), mp.mpf(1), mp.mpf('2.2')], [1, 1, 1], [0, 0, 0])
    out['hexA_w123_power'] = run('hexA_w123_power', [mp.mpf(0), -mp.pi/3, -2*mp.pi/3], [1, 2, 3],
                                 [0, mp.log(2), mp.log(3)])
    json.dump(out, open(f's06_out_{mp.mp.dps}.json', 'w'), indent=1)
