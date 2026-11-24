"""Reproduce CX2's exact checks and independent polar finite differences.

Run with the project's Python: review.py exact | formula | basis | numeric.
All inputs outside this directory are opened read-only; no K1e code is imported.
"""
from pathlib import Path
from math import comb
import argparse
import pickle
import re
import time

import mpmath as mp
import sympy as sp


SOURCE = Path(__file__).resolve().parent.parent / "K1e_hexagon_closed"
VERTICES = [sp.Matrix(v) for v in [(2, 0), (3, 2), (0, 1)]]
WEIGHTS = list(map(sp.Integer, [2, 5, 7]))
ORDER = ["Li3(x)", "Li3(x/(x-1))", "Li2(x)*ln(1-x)",
         "ln^3(1-x)", "Li2(x)", "ln^2(1-x)", "ln(1-x)",
         "zeta3", "zeta2*ln(1-x)", "zeta2", "1"]


def cross(a, b):
    return a[0]*b[1] - a[1]*b[0]


def number(x):
    x = sp.Rational(x)
    return mp.mpf(int(x.p))/int(x.q)


def exact():
    c, v = [1], [1]
    for n in range(1, 62):
        nc = (2**n if n % 2 == 0 else 0) - sum(
            comb(n, k)*c[k]*c[n-k] for k in range(1, n))
        nv = 2**(2*n-1) - sum(
            comb(n, k)*v[k]*v[n-k] for k in range(1, n))
        assert nc % 2 == nv % 2 == 0
        c.append(nc//2)
        v.append(nv//2)
        assert c[n] % 2 == 0 and v[n] % 2 == 1
        assert v[n] == sum(comb(n, k)*c[k] for k in range(n+1))
    print("n,V_n,V_n mod 2")
    for n, value in enumerate(v):
        print(n, value, value % 2, sep=",")
    print("c_0 through c_8:", c[:9])
    print("Catalan coefficients a_1 through a_8:",
          [(-1)**(j-1)*2*sp.catalan(j-1) for j in range(1, 9)])
    for n, expected in [(3, sp.Integer(1)), (5, -sp.Rational(3, 4))]:
        omega = sp.Rational(4*(-1)**n, n)*(1-sp.Rational(v[n], 2**(n-1)))
        print("omega_" + str(n), "=", omega)
        assert omega == expected
    print("ALL EXACT CHECKS PASS through n=61; n=1 gives omega_1=0.")


def geometry(verbose=True):
    matrix = sp.Matrix([[abs(cross(a, b)) for b in VERTICES] for a in VERTICES])
    effective = list(matrix.inv()*sp.ones(3, 1))
    k = [effective[j]/WEIGHTS[j] for j in range(3)]
    verts = VERTICES + [-v for v in VERTICES]
    turns = [cross(verts[(j+1) % 6]-verts[j],
                   verts[(j+2) % 6]-verts[(j+1) % 6]) for j in range(6)]
    assert all(c > 0 for c in effective) and all(t > 0 for t in turns)
    assert matrix*sp.Matrix(effective) == sp.ones(3, 1)
    edges = []
    for q in range(3):
        a, b = verts[q:q+2]
        d = b-a
        j0 = next(j for j, v in enumerate(VERTICES) if cross(v, a) == 0)
        j1 = next(j for j, v in enumerate(VERTICES) if cross(v, b) == 0)
        jo = next(j for j in range(3) if j not in (j0, j1))
        e = dict(A0=effective[j0]*abs(cross(VERTICES[j0], b)),
                 A1=effective[j1]*abs(cross(VERTICES[j1], a)),
                 B=effective[jo]*abs(cross(VERTICES[jo], a)),
                 r=cross(VERTICES[jo], b)/cross(VERTICES[jo], a),
                 L2=d.dot(d), q0=a.dot(d), V02=a.dot(a))
        e["kap"] = e["r"]-1
        assert e["A0"] == 1-e["B"]*e["r"] > 0
        assert e["A1"]+e["B"] == 1 and e["r"] > 0
        e["lA0"] = sp.log(e["A0"]/WEIGHTS[j0])
        e["lA1"] = sp.log(e["A1"]/WEIGHTS[j1])
        e["lB"] = sp.log(e["B"]/WEIGHTS[jo])
        e["Z"] = sp.simplify(6*e["A0"]*e["A1"]*
            (e["L2"]*e["V02"]-e["q0"]**2)*
            (3*e["A1"]*e["L2"]+3*e["B"]*e["kap"]*e["q0"]-2*e["L2"])/
            e["L2"]**sp.Rational(5, 2))
        edges.append(e)
    if verbose:
        print("vertices:", [list(v) for v in VERTICES], "+ negatives")
        print("weights:", WEIGHTS, "effective weights w*k:", effective, "k:", k)
        print("N(V_q):", list(matrix*sp.Matrix(effective)), "convexity turns:", turns)
        for q, e in enumerate(edges):
            print("edge", q+1, {key: str(val) for key, val in e.items()})
    return verts, k, edges


def coefficients():
    # This is the supplied local symbolic artifact, not executable K1e code.
    data = pickle.loads((SOURCE / "Earc_basis_clean.pkl").read_bytes())
    assert set(data) == set(ORDER)
    symbols = sp.symbols("A1 B kap q0 V02 lA0 lA1 lB L2")
    expressions = [data[name][0]/data[name][1] for name in ORDER]
    return data, sp.lambdify(symbols, expressions, "mpmath")


def atoms(r, sign):
    x, ell = 1-r, sign*mp.log(r)
    li2 = mp.polylog(2, x)
    return [mp.polylog(3, x), mp.polylog(3, x/(x-1)), li2*ell,
            ell**3, li2, ell**2, ell, mp.zeta(3), mp.zeta(2)*ell,
            mp.zeta(2), mp.mpf(1)]


def formula():
    _, _, edges = geometry()
    data, evaluate = coefficients()
    totals = {1: mp.mpf(0), -1: mp.mpf(0)}
    A1, B, kap, r, L2 = sp.symbols("A1 B kap r L2")
    blocks = re.split(r"^===== coefficient of (.*?)  \([^\n]*=====\n",
                      (SOURCE / "formula_table.txt").read_text(), flags=re.M)
    assert len(blocks) == 1+2*len(ORDER)
    for i in range(1, len(blocks), 2):
        name, body = blocks[i:i+2]
        table = sum(sp.sympify(m)*sp.sympify(c) for m, c in
                    re.findall(r"^  \[(.*?)\] \* \( (.*?) \)$", body, re.M))
        num, den = data[name]
        assert den == L2**sp.Rational(5, 2)
        assert sp.expand(table-num.subs({A1: 1-B, kap: r-1})) == 0
    print("Human table vs pickle: all 11 numerator differences exactly zero.")
    for q, e in enumerate(edges):
        args = [number(e[s]) if not s.startswith("l") else
                mp.log(number(sp.exp(e[s])))
                for s in ["A1", "B", "kap", "q0", "V02", "lA0", "lA1", "lB", "L2"]]
        co = evaluate(*args)
        zexact = sp.simplify((data["zeta3"][0]/data["zeta3"][1]).subs(
            {sp.Symbol(s): e[s] for s in ["A1", "B", "kap", "q0", "V02", "L2"]}))
        assert sp.simplify(zexact-e["Z"]) == 0
        print("edge", q+1, "r=", e["r"], "exact zeta3 coefficient=Z_e=", zexact)
        for sign in (1, -1):
            val = mp.fsum(c*a for c, a in zip(co, atoms(number(e["r"]), sign)))
            totals[sign] += 2*val
            print("  ell=" + ("+" if sign == 1 else "-") + "log(r):", mp.nstr(val, 90))
    for sign, total in totals.items():
        print("FORMULA TOTAL ell=" + ("+" if sign == 1 else "-") + "log(r):", mp.nstr(total, 90))


def basis():
    # Sampling cannot establish independence over Q(r); see RESULT.md.
    points = [sp.Rational(x) for x in ["1/4", "1/3", "2/5", "1/2", "2/3", "3/4",
                                      "3/2", "2", "3", "4", "5", "6"]]
    full = mp.matrix([atoms(number(r), 1) for r in points])
    residual = max(abs(full[i, 7]-mp.zeta(3)*full[i, 10]) for i in range(len(points)))
    print("sample r:", points)
    print("zeta3 column - zeta3*one column max residual:", mp.nstr(residual, 12))
    _, singular, _ = mp.svd(full)
    print("11-column singular values:", [mp.nstr(x, 12) for x in singular])
    print("weight-3 word coefficient matrix, columns baa,bab,bba,bbb:")
    words = sp.Matrix([[-1, 0, 0, 0], [1, -1, -1, 1],
                       [0, -1, -2, 0], [0, 0, 0, 6]])
    print(words)
    print("word matrix determinant:", words.det())
    keep = [0, 1, 2, 3, 4, 5, 6, 10]
    for degree in [0, 1, 2]:
        n = len(keep)*(degree+1)
        rows = []
        for i in range(1, n+1):
            r = mp.mpf(i)/(n+1-i)
            f = atoms(r, 1)
            rows.append([r**j*f[k] for j in range(degree+1) for k in keep])
        m = mp.matrix(rows)
        # Scaling avoids mistaking large column ranges for rank loss.
        for j in range(n):
            scale = mp.sqrt(mp.fsum(m[i, j]**2 for i in range(n)))
            for i in range(n):
                m[i, j] /= scale
        _, sv, _ = mp.svd(m)
        print("8-function polynomial-coefficient test degree", degree,
              "matrix size", n, "normalized determinant", mp.nstr(mp.det(m), 20),
              "s_min/s_max", mp.nstr(sv[n-1]/sv[0], 20))


def numeric(h0, levels):
    verts, k, _ = geometry()
    started = time.monotonic()
    perimeters = {}
    arc_values = {}
    arc_integrators = []
    for q in range(3):
        a, b = verts[q:q+2]
        a, b = [[number(x) for x in v] for v in (a, b)]
        area = cross(a, b)
        ends = [(number(k[j])*cross([number(x) for x in v], a),
                 number(k[j])*cross([number(x) for x in v], b),
                 number(k[j])*sum(number(v[z])*a[z] for z in range(2)),
                 number(k[j])*sum(number(v[z])*b[z] for z in range(2)),
                 number(WEIGHTS[j])) for j, v in enumerate(VERTICES)]
        cache = {}

        def integral(p, a=a, b=b, area=area, ends=ends, cache=cache):
            def integrand(t):
                if t not in cache:
                    u = 1-t
                    y = [u*a[j]+t*b[j] for j in range(2)]
                    radius2 = sum(z*z for z in y)
                    radius = mp.sqrt(radius2)
                    terms = []
                    for g0, g1, gp0, gp1, w in ends:
                        g = u*g0+t*g1
                        s = abs(g)/radius
                        # Exact endpoints have measure zero. Tanh-sinh skips them.
                        assert s > 0
                        terms.append((mp.log(s), w*s,
                                      w*mp.sign(g)*(u*gp0+t*gp1)/radius))
                    cache[t] = (area/radius2, terms)
                jacobian, terms = cache[t]
                powers = [mp.exp((p-1)*slog) for slog, _, _ in terms]
                H = mp.fsum(a*e for (_, a, _), e in zip(terms, powers))
                T = mp.fsum(b*e for (_, _, b), e in zip(terms, powers))
                # H_theta = p*T; rho_theta/rho = -T/H.
                return mp.exp(-mp.log(H)/p)*mp.sqrt(1+(T/H)**2)*jacobian
            return mp.quad(integrand, [0, mp.mpf("0.5"), 1], method="tanh-sinh")
        arc_integrators.append(integral)

    def perimeter(delta):
        if delta not in perimeters:
            arcs = [f(1+delta) for f in arc_integrators]
            val = 2*mp.fsum(arcs)
            perimeters[delta] = val
            arc_values[delta] = mp.matrix(arcs)
            print("P(1" + ("+" if delta >= 0 else "") + mp.nstr(delta, 14) + ") =",
                  mp.nstr(val, 85), "arcs=", [mp.nstr(v, 85) for v in arcs], flush=True)
        return perimeters[delta]

    baseline = perimeter(mp.mpf(0))
    expected = 2*(2*mp.sqrt(5)+mp.sqrt(10))
    print("P(1) exact polygon check |diff|:", mp.nstr(abs(baseline-expected), 12), flush=True)
    assert abs(baseline-expected) < mp.mpf("1e-60")
    previous = []
    previous_arcs = []
    for level in range(levels):
        h = h0/2**level
        d3 = (perimeter(2*h)-2*perimeter(h)+2*perimeter(-h)-perimeter(-2*h))/(2*h**3)
        row = [d3]
        arc_row = [(arc_values[2*h]-2*arc_values[h]+2*arc_values[-h]-arc_values[-2*h])/(2*h**3)]
        for j in range(1, level+1):
            row.append((4**j*row[j-1]-previous[j-1])/(4**j-1))
            arc_row.append((4**j*arc_row[j-1]-previous_arcs[j-1])/(4**j-1))
        print("RICHARDSON level", level, "h=", mp.nstr(h, 12),
              "row=", [mp.nstr(x, 75) for x in row], flush=True)
        if previous:
            print("diagonal change:", mp.nstr(abs(row[-1]-previous[-1]), 12), flush=True)
        previous = row
        previous_arcs = arc_row
    print("DIRECT POLAR FD FINAL:", mp.nstr(previous[-1], 85), flush=True)
    print("DIRECT POLAR FD ARCS:", [mp.nstr(v, 85) for v in previous_arcs[-1]], flush=True)
    _, evaluate = coefficients()
    _, _, edges = geometry(False)
    for sign in (1, -1):
        total = mp.mpf(0)
        for q, e in enumerate(edges):
            args = [number(e[s]) if not s.startswith("l") else mp.log(number(sp.exp(e[s])))
                    for s in ["A1", "B", "kap", "q0", "V02", "lA0", "lA1", "lB", "L2"]]
            co = evaluate(*args)
            arc = mp.fsum(c*a for c, a in zip(co, atoms(number(e["r"]), sign)))
            total += 2*arc
            if sign == -1:
                print("corrected formula-FD edge", q+1, ":",
                      mp.nstr(arc-previous_arcs[-1][q], 20), flush=True)
        print("formula ell=" + ("+" if sign == 1 else "-") + "log(r) total:",
              mp.nstr(total, 85), "formula-FD:", mp.nstr(total-previous[-1], 20), flush=True)
        if sign == -1:
            assert abs(total-previous[-1]) < mp.mpf("1e-20")
    print("elapsed seconds:", round(time.monotonic()-started, 2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["exact", "formula", "basis", "numeric"])
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--h", default="0.001")
    parser.add_argument("--levels", type=int, default=6)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    if args.mode == "numeric":
        assert mp.mpf(args.h) > 0 and 2*mp.mpf(args.h) < mp.mpf("0.5")
        assert args.levels >= 2 and args.dps >= 75
        numeric(mp.mpf(args.h), args.levels)
    else:
        globals()[args.mode]()
