"""Second algorithm (mpmath.pslq) on the 16-atom support of the c_1 relation, at two precisions (data dps 600 and dps 560).
Atoms: lam^4, lam^3 L, lam^2 L^2, lam L^3, L^4, pi^2 lam^2, pi^2 lam L, pi^2 L^2, pi^4, z3 lam, z3 L, Li3(b) lam, Li3(b) L, Li3(-b) lam, Li3(-b) L, Li3(1/sqrt2) L."""
from id_lib import *
names = ["lam^4", "lam^3 L", "lam^2 L^2", "lam L^3", "L^4", "pi2 lam^2", "pi2 lam L", "pi2 L^2", "pi4", "z3 lam", "z3 L",
         "Li3b lam", "Li3b L", "Li3mb lam", "Li3mb L", "Li3h L"]
for dps, tol in ((600, 570), (560, 530)):
    set_dps(dps)
    I = load_I(dps); c1 = combo(I, "c1"); A = atoms()
    rel = mp.pslq([c1] + [A[n] for n in names], maxcoeff=10**6, maxsteps=10**7, tol=mp.mpf(10)**(-tol))
    res = sum(r * x for r, x in zip(rel, [c1] + [A[n] for n in names])) if rel else None
    print(f"dps {mp.mp.dps} tol 1e-{tol}: relation {rel}")
    print("   residual", mp.nstr(res, 5), " (mp.mp.dps", mp.mp.dps, ")")
