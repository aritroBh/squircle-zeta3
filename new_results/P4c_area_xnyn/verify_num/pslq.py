# Identify the freely fitted constants (config A, rho=1) by integer relation search, with gamma_E and pi allowed in the basis.
import json, mpmath as mp
res = json.load(open('fit_rho1.json'))
fA = res['A']['fit']
for dps, keys in [(45, ['(1, 1)', '(1, 0)', '(2, 2)', '(2, 1)', '(2, 0)'])]:
    mp.mp.dps = dps
    lam, z2, z3, g, s2 = mp.log(2), mp.pi**2/6, mp.zeta(3), mp.euler, mp.sqrt(2)
    bas = [('1', mp.mpf(1)), ('sqrt2', s2), ('lam', lam), ('lam^2', lam**2), ('lam^3', lam**3), ('z2', z2), ('z2*lam', z2*lam), ('z3', z3), ('gammaE', g), ('pi', mp.pi)]
    names = {'(1, 1)': 'beta', '(1, 0)': 'gamma', '(2, 2)': 'b3', '(2, 1)': 'c3', '(2, 0)': 'd3'}
    for k in keys:
        x = mp.mpf(fA[k])
        rel = mp.pslq([x] + [v for _, v in bas], maxcoeff=10**4, maxsteps=10**6)
        if rel is None: print(names[k], 'no relation'); continue
        c0 = rel[0]
        terms = [f"{mp.mpf(-c)/c0}*{nm}" for c, (nm, _) in zip(rel[1:], bas) if c != 0]
        print(f"{names[k]} = " + ' + '.join(terms), '   [relation vector', rel, ']')
