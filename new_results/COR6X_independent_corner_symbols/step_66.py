t = r'''
import sympy as sp, mpmath as mp
mp.mp.dps = 50
xx = sp.Symbol('xx')
for n in [7, 9, 11, 3, 5]:
    c0 = sp.sin(sp.pi/n)
    mp0 = sp.minpoly(c0, xx)
    # true field: Q(sin(pi/n)) has degree phi(2n)/2 for n != ... compute via cos(2pi/n)
    c2 = sp.cos(2*sp.pi/n)
    mp2 = sp.minpoly(c2, xx)
    # check sin(pi/n) satisfies mp2 via sin^2 = (1-cos2)/2
    v = mp.nstr(mp.polyval([float(c) for c in sp.Poly(mp2, xx).all_coeffs()], 1 - 2*mp.sin(mp.pi/n)**2), 3) if False else None
    # numerical check: sin^2(pi/n) = (1-cos(2pi/n))/2
    lhs = mp.sin(mp.pi/n)**2
    rhs = (1 - mp.cos(2*mp.pi/n))/2
    print(f"n={n}: deg minpoly(sin)={sp.degree(mp0)}, deg minpoly(cos2)={sp.degree(mp2)}, sin^2==...: {mp.nstr(abs(lhs-rhs),3)}")
    print("   minpoly(sin):", mp0.as_expr())
    print("   minpoly(cos(2pi/n)):", mp2.as_expr())
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t], capture_output=True, text=True)
print(r.stdout); print("ERR:", (r.stderr or "")[-500:])
