"""
H1: large-order behaviour of the cross-polytope-end coefficients (numerical).

Data: hi_coeffs_n{n}_K80.json (05_high_order_numeric.py; eta- and d-coefficients of S_n/S_n(1), 150 digits;
orders <= 12 agree with the exact Routes A/B).
Prediction (heuristic, result.md Sec. 6): the factorial divergence comes only from the bending factor B_n, and the
nearest Borel singularities sit where the integrand sqrt(mean_i exp(2 eta m_i)) first vanishes along the
'edge' direction (one coordinate -> 0, which has the slowest exponential tail):
     t_edge = ( -log(n-1) +- i pi ) / 2 ,     |t_edge| = sqrt(pi^2 + log^2(n-1)) / 2.
For n = 2 this is +- i pi/2, the known diamond value (D11, R2).  The 'vertex' direction (two coordinates -> 0)
gives t_vertex = log 2 +- i pi for n = 3 (farther away).
Tests: (a) Pade poles of the Borel transform sum (a_k/k!) t^k; (b) root test |a_k/k!|^{-1/k}.
"""
import json, sys
import mpmath as mp

mp.mp.dps = 140
for n in (2, 3, 4, 5):
    H = json.load(open(f'hi_coeffs_n{n}_K80.json'))
    for var in ('eta', 'd'):
        a = [mp.mpf(x) for x in H[var]]
        K = len(a) - 1
        b = [a[k] / mp.factorial(k) for k in range(K + 1)]
        te = (-mp.log(n - 1) + 1j * mp.pi) / 2
        print(f'===== n = {n}, variable {var}: predicted edge singularity t = {mp.nstr(te, 12)}, |t| = {mp.nstr(abs(te), 12)}')
        for (L, M) in ((30, 30), (38, 38), (40, 40)):
            p, q = mp.pade(b[:L + M + 1], L, M)
            roots = mp.polyroots(q[::-1], maxsteps=400, extraprec=400)
            roots = sorted(roots, key=lambda r: abs(r))
            # keep poles that are not cancelled by a nearby zero of p
            proots = mp.polyroots(p[::-1], maxsteps=400, extraprec=400) if L > 0 else []
            good = []
            for r in roots:
                dz = min([abs(r - s) for s in proots]) if proots else mp.inf
                if dz > 1e-8 * max(1, abs(r)):
                    good.append(r)
            print(f'  Pade [{L}/{M}] nearest genuine Borel poles: ' +
                  ', '.join(mp.nstr(r, 10) for r in good[:6]))
        rt = [mp.nstr(abs(b[k]) ** (-mp.mpf(1) / k), 8) for k in (40, 50, 60, 70, 80)]
        print(f'  root test |a_k/k!|^(-1/k) at k = 40,50,60,70,80: {rt}')
