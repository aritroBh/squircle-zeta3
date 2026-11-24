"""
H1: identification of c_1, c_2, c_3 (coefficients of S(1+d)/S(1), n = 3) WITHOUT the Gamma-function
algebra of Routes A/B.

Start from the simplex integral (checked against the direct graph quadrature N1 to 40 digits, file
direct_S_dps40.json):  S(p) = 8 eps^2 Int_Delta (uvw)^{-eta} sqrt(u^{2eta}+v^{2eta}+w^{2eta}) du dv.
Expand the integrand pointwise in eta: (uvw)^{-eta} sqrt(sum u_i^{2eta}) = sqrt3 * sum_k eta^k p_k(log u, log v, log w),
p_k polynomials (sympy), integrate each p_k over the simplex by 2-D tanh-sinh at 70 digits,
then convert eta -> d and run PSLQ on small bases.
Basis choice: by the structure theorem the d^k coefficient has weight <= k in single zetas, so
{1, pi^2, zeta(3)} is the natural basis for k <= 3.  To test for 'edge/vertex' constants we also run
an enlarged basis containing log2, log3, pi*sqrt3, sqrt3, L(chi_-3, 2) (= Cl2(2pi/3)*2/sqrt3... see code),
Catalan G, pi*log2: a genuine relation must not use them.
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 70
eta = sp.Symbol('eta')
L1, L2, L3 = sp.symbols('L1 L2 L3')
KMAX = 4
expr = sp.exp(-eta * (L1 + L2 + L3)) * sp.sqrt((sp.exp(2 * eta * L1) + sp.exp(2 * eta * L2) + sp.exp(2 * eta * L3)) / 3)
ser = sp.series(expr, eta, 0, KMAX + 1).removeO()
polys = [sp.expand(ser.coeff(eta, k)) for k in range(KMAX + 1)]
for k, pk in enumerate(polys):
    print(f'p_{k} = {pk}')
fns = [sp.lambdify((L1, L2, L3), pk, 'mpmath') for pk in polys]


def simplex_int(fn):
    def inner(u):
        lu = mp.log(u)
        return (1 - u) * mp.quad(lambda t: fn(lu, mp.log((1 - u) * t), mp.log((1 - u) * (1 - t))), [0, 1])
    return mp.quad(inner, [0, 1])


J = []
for k, fn in enumerate(fns):
    val = simplex_int(fn)
    J.append(val)
    print(f'J_{k} = Int_Delta p_{k} = {mp.nstr(val, 60)}')
# S/S(1) = 8 (1-eta)^2 sqrt3 * sum J_k eta^k / (4 sqrt3) = 2 (1-eta)^2 sum J_k eta^k
ce = [2 * J[k] - 4 * (J[k - 1] if k >= 1 else 0) + 2 * (J[k - 2] if k >= 2 else 0) for k in range(KMAX + 1)]
# eta = d/(1+d): coefficients in d
dsym = sp.Symbol('d')
cd = [mp.mpf(0)] * (KMAX + 1)
etad = sp.series(dsym / (1 + dsym), dsym, 0, KMAX + 1).removeO()
for k in range(KMAX + 1):
    pw = sp.Poly(sp.expand(etad ** k), dsym)
    for (e,), coef in zip(pw.monoms(), pw.coeffs()):
        if e <= KMAX:
            cd[e] += ce[k] * mp.mpf(sp.Rational(coef).p) / mp.mpf(sp.Rational(coef).q)
for k in range(KMAX + 1):
    print(f'c_{k} (d-coefficient of S/S(1)) = {mp.nstr(cd[k], 60)}')

z3 = mp.zeta(3)
Lchi = (mp.zeta(2, mp.mpf(1) / 3) - mp.zeta(2, mp.mpf(2) / 3)) / 9  # L(chi_{-3}, 2)
print('\n--- PSLQ, natural basis ---')
for k, basis_names, basis in [
    (1, ['c1', '1'], [cd[1], 1]),
    (2, ['c2', '1', 'pi^2'], [cd[2], 1, mp.pi ** 2]),
    (3, ['c3', '1', 'pi^2', 'zeta3'], [cd[3], 1, mp.pi ** 2, z3]),
    (4, ['c4', '1', 'pi^2', 'zeta3', 'pi^4'], [cd[4], 1, mp.pi ** 2, z3, mp.pi ** 4]),
]:
    for dig in (50, 60):
        mp.mp.dps = dig
        rel = mp.pslq([mp.mpf(b) for b in basis], maxcoeff=10 ** 6, maxsteps=10 ** 6)
        res = sum(r * b for r, b in zip(rel, basis)) if rel else None
        print(f'  c_{k} @ {dig} digits: relation {dict(zip(basis_names, rel)) if rel else None}, '
              f'residual {mp.nstr(res, 3) if rel else "-"}, max|coef| {max(abs(r) for r in rel) if rel else "-"}')
    mp.mp.dps = 70
print('\n--- PSLQ, enlarged basis (edge/vertex-type constants allowed) ---')
names = ['c', '1', 'pi^2', 'zeta3', 'log2', 'log3', 'sqrt3', 'pi*sqrt3', 'sqrt3*L(chi-3,2)', 'Catalan', 'pi*log2', 'pi^2*sqrt3']
for k in (2, 3):
    basis = [cd[k], 1, mp.pi ** 2, z3, mp.log(2), mp.log(3), mp.sqrt(3), mp.pi * mp.sqrt(3), mp.sqrt(3) * Lchi, mp.catalan,
             mp.pi * mp.log(2), mp.pi ** 2 * mp.sqrt(3)]
    mp.mp.dps = 68
    rel = mp.pslq([mp.mpf(b) for b in basis], maxcoeff=10 ** 5, maxsteps=10 ** 7)
    res = sum(r * b for r, b in zip(rel, basis)) if rel else None
    print(f'  c_{k} @ 68 digits: {dict((nm, r) for nm, r in zip(names, rel) if r != 0) if rel else None}, '
          f'residual {mp.nstr(res, 3) if rel else "-"}')
mp.mp.dps = 70
print('\nexact predictions: c1 = 1, c2 = 1 - pi^2/9, c3 = 1 + pi^2/9 - 16 zeta3/9, c4 = 1 - 2pi^2/9 + 32 zeta3/9 - 16 pi^4/405')
pred = [1, 1, 1 - mp.pi ** 2 / 9, 1 + mp.pi ** 2 / 9 - 16 * z3 / 9, 1 - 2 * mp.pi ** 2 / 9 + 32 * z3 / 9 - 16 * mp.pi ** 4 / 405]
for k in range(KMAX + 1):
    print(f'  c_{k}: numeric - exact = {mp.nstr(cd[k] - pred[k], 3)}')
