# s13: the decomposition for ANY number of normals (m >= 4: several non-local normals per edge, products
# N_i N_j appear).  P'''(1) = Z zeta(3) + Q0 + Q2 zeta(2) + 6*NL, with NL a convergent multiple power series
# whose terms use only Beta moments with a+b <= 2 (rationals and zeta(2)).  Checked against D15b's polar total.
import sys, json, pickle, mpmath as mp, sympy as sp
sys.path.insert(0, 'research/08_theory/D15b_polygon_law')
import poly_coeffs as pc
from edge_num import edges, zeta3_edge_law
from s03_beta_table import J, z2, z3

mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
D = pickle.load(open('s01_E.pkl', 'rb')); E3 = D['E'][3]
t, al, be, L, h, s0 = sp.symbols('t alpha beta L h s0', real=True)
l0, l1 = sp.symbols('l0 l1'); tau = sp.symbols('tau1:4'); gam = sp.symbols('gam1:4')
X, Y, mu0, mu1 = sp.symbols('X Y mu0 mu1')
HC = {1: [mp.mpf(0)], 2: [mp.mpf(0)]}
def Hn(q, r):
    arr = HC[r]
    while len(arr) <= q: arr.append(arr[-1] + mp.mpf(1)/len(arr)**r)
    return arr[q]
def Jnum2(n, k, a, b):
    assert a + b <= 2
    B = mp.factorial(n)*mp.factorial(k)/mp.factorial(n+k+1)
    fx = Hn(n, 1) - Hn(n+k+1, 1); fy = Hn(k, 1) - Hn(n+k+1, 1)
    fxx = -Hn(n, 2) + Hn(n+k+1, 2); fyy = -Hn(k, 2) + Hn(n+k+1, 2); fxy = -(mp.zeta(2) - Hn(n+k+1, 2))
    return {(0, 0): B, (1, 0): B*fx, (0, 1): B*fy, (2, 0): B*(fx**2 + fxx), (0, 2): B*(fy**2 + fyy),
            (1, 1): B*(fx*fy + fxy)}[(a, b)]

def build_poly(q):
    A = sp.symbols(f'A1:{q+1}'); B = sp.symbols(f'B1:{q+1}'); MU = sp.symbols(f'mu2:{q+2}'); NN = sp.symbols(f'N1:{q+1}')
    sub = {l0: X + mu0, l1: Y + mu1}
    for k in range(1, 4):
        sub[tau[k-1]] = sum((A[i] + B[i]*t)*(MU[i] + NN[i])**k for i in range(q))
        sub[gam[k-1]] = sum(B[i]*(MU[i] + NN[i])**k for i in range(q))
    P = sp.Poly(sp.expand(E3.xreplace(sub)), t, X, Y, *NN)
    syms = (al, be, L, h, s0, mu0, mu1) + A + B + MU
    return [(mon, sp.lambdify(syms, c, 'mpmath')) for mon, c in P.terms()]

def nseries(r, nmax):
    """power series of N = log(1+(r-1)t): returns (mode, coefficient list) in t or in s = 1-t"""
    if abs(r - 1) <= abs(1/r - 1):
        x = r - 1; return 't', [mp.mpf(0)] + [(-1)**(q+1)*x**q/q for q in range(1, nmax+1)]
    xp = 1/r - 1; return 's', [mp.log(r)] + [(-1)**(q+1)*xp**q/q for q in range(1, nmax+1)]
def spow(c, j, nmax):
    out = [mp.mpf(1)] + [mp.mpf(0)]*nmax
    for _ in range(j):
        out = [mp.fsum([out[i]*c[q-i] for i in range(q+1)]) for q in range(nmax+1)]
    return out

def run(name, bs, ws, lams):
    m = len(bs); q = m - 2
    terms = build_poly(q)
    Z, parts = zeta3_edge_law(bs, ws)
    Q0 = Q2 = Q3 = NL = NLq = mp.mpf(0); worst_edge = mp.mpf(0)
    for e in edges(bs, ws, lams):
        j0, j1 = e['j0'], e['j1']; nl = [i for i in range(m) if i not in (j0, j1)]
        Ai, Bi, MUi, ri = [], [], [], []
        for i in nl:
            sg = mp.sign(e['c0'][i]); Ai.append(ws[i]*abs(e['c0'][i])); Bi.append(ws[i]*sg*(e['c1'][i] - e['c0'][i]))
            MUi.append(mp.log(abs(e['c0'][i])) + lams[i]); ri.append(e['c1'][i]/e['c0'][i])
        m0 = mp.log(abs(e['c1'][j0])) + lams[j0]; m1 = mp.log(abs(e['c0'][j1])) + lams[j1]
        vals = (e['alpha'], e['beta'], e['L'], e['h'], e['s0'], m0, m1) + tuple(Ai) + tuple(Bi) + tuple(MUi)
        coef = [(mon, f(*vals)) for mon, f in terms]
        rho = max(min(abs(r - 1), abs(1/r - 1)) for r in ri)
        nmax = int((mp.mp.dps + 8)/(-mp.log10(rho))) + 25
        ser = [nseries(r, nmax) for r in ri]
        lq = [mp.mpf(0)]*3; nls = mp.mpf(0)
        for mon, c in coef:
            n, a, b = mon[:3]; js = mon[3:]
            if sum(js) == 0:
                ex = J(n, 0, a, b)
                lq[0] += c*mp.mpf(str(sp.N(ex.subs({z2: 0, z3: 0}), mp.mp.dps + 10)))
                lq[1] += c*mp.mpf(str(sp.N(ex.coeff(z2), mp.mp.dps + 10)))
                lq[2] += c*mp.mpf(str(sp.N(ex.coeff(z3), mp.mp.dps + 10)))
                continue
            # product of series in t and series in s
            tser = [mp.mpf(1)] + [mp.mpf(0)]*nmax; sser = [mp.mpf(1)] + [mp.mpf(0)]*nmax
            for (mode, cs), j in zip(ser, js):
                if j == 0: continue
                pw = spow(cs, j, nmax)
                if mode == 't': tser = [mp.fsum([tser[i]*pw[k-i] for i in range(k+1)]) for k in range(nmax+1)]
                else: sser = [mp.fsum([sser[i]*pw[k-i] for i in range(k+1)]) for k in range(nmax+1)]
            acc = mp.mpf(0)
            for qa in range(nmax+1):
                if tser[qa] == 0: continue
                for qb in range(nmax+1 - qa):
                    if sser[qb] == 0: continue
                    acc += tser[qa]*sser[qb]*Jnum2(n + qa, qb, a, b)
            nls += c*acc
        Nf = [lambda tt, r=r: mp.log(1 + (r - 1)*tt) for r in ri]
        nlq = mp.quad(lambda tt: mp.fsum([c*tt**mon[0]*mp.log(tt)**mon[1]*mp.log(1-tt)**mon[2]*
                      mp.fprod([Nf[i](tt)**mon[3+i] for i in range(q)]) for mon, c in coef if sum(mon[3:]) > 0]),
                      [0, mp.mpf(1)/2, 1])
        worst_edge = max(worst_edge, abs(nls - nlq))
        Q0 += 6*lq[0]; Q2 += 6*lq[1]; Q3 += 6*lq[2]; NL += 6*nls
    P3_dec = Q0 + Q2*mp.zeta(2) + Q3*mp.zeta(3) + NL
    P3_polar = 6*pc.P_coeffs_general(bs, ws, lams, kmax=3)[3]
    print(f"{name} (m={m}): Z edge law = {mp.nstr(Z, 25)}; 6*local z3-coef = {mp.nstr(Q3, 25)}; "
          f"max edge |NL series - NL quad| = {mp.nstr(worst_edge, 3)}")
    print(f"   P'''(1): decomposed {mp.nstr(P3_dec, 30)}  polar(D15b) {mp.nstr(P3_polar, 30)}  |diff| {mp.nstr(abs(P3_dec - P3_polar), 3)}")
    print(f"   remainder: Q0 = {mp.nstr(Q0, 20)}, Q2 = {mp.nstr(Q2, 20)}, 6NL = {mp.nstr(NL, 20)}")
    return dict(Z=str(Z), Q3=str(Q3), Q0=str(Q0), Q2=str(Q2), NL=str(NL), dec=str(P3_dec), polar=str(P3_polar))

if __name__ == '__main__':
    out = {}
    out['reg_m4'] = run('regular octagon', [-(i*mp.pi/4) for i in range(4)], [1]*4, [0]*4)
    out['oct_generic'] = run('generic octagon', [mp.mpf(0), mp.mpf('0.7'), mp.mpf('1.6'), mp.mpf('2.5')],
                             [1, mp.mpf('1.5'), mp.mpf('0.8'), mp.mpf('1.2')], [0, 0, 0, 0])
    json.dump(out, open(f's13_out_{mp.mp.dps}.json', 'w'), indent=1)
