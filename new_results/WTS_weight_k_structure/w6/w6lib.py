"""shared helpers for the w6 identification scripts (imports only json/os/mpmath; precision is always set by the caller AFTER import)."""
import os, json
import mpmath as mp
W = os.path.dirname(os.path.abspath(__file__))

def coef(run, j, last=False):
    """c_j of W_6(k) from the exact interpolation of the high-precision run `run` (dps of that run), as a decimal string."""
    d = json.load(open(os.path.join(W, f"coef_{run}.json")))
    return d["coefficients_last_fit" if last else "coefficients"][j]

def consts():
    """alternating-MZV style constants at the CURRENT mp.mp.dps (caller must have set it; S6 cache has ~275 digits)."""
    c = json.load(open(os.path.join(W, "consts_270.json")))
    lam = mp.log(2); z3 = mp.zeta(3); z5 = mp.zeta(5); pi2 = mp.pi ** 2
    h = mp.mpf(1) / 2
    return dict(lam=lam, z3=z3, z5=z5, pi2=pi2, Li4=mp.polylog(4, h), Li5=mp.polylog(5, h), Li6=mp.polylog(6, h),
                S6=mp.mpf(c["S6"]))

def basis13(C):
    lam, z3, z5, pi2, Li4, Li5, Li6, S6 = (C[x] for x in ("lam", "z3", "z5", "pi2", "Li4", "Li5", "Li6", "S6"))
    return [("z3^2", z3 ** 2), ("z5 lam", z5 * lam), ("pi^6", pi2 ** 3), ("pi^2 z3 lam", pi2 * z3 * lam),
            ("z3 lam^3", z3 * lam ** 3), ("pi^4 lam^2", pi2 ** 2 * lam ** 2), ("pi^2 lam^4", pi2 * lam ** 4),
            ("lam^6", lam ** 6), ("pi^2 Li4(1/2)", pi2 * Li4), ("lam^2 Li4(1/2)", lam ** 2 * Li4),
            ("lam Li5(1/2)", lam * Li5), ("Li6(1/2)", Li6), ("S6=zeta(-5,-1)", S6)]

def pslq_try(target, basis, D, maxcoeff=10 ** 12, maxsteps=10 ** 6):
    """PSLQ on [target]+basis at the current precision; returns (relation or None, residual, tol)."""
    tol = mp.mpf(10) ** (-(D - 25))
    vec = [target] + [v for _, v in basis]
    rel = mp.pslq(vec, tol=tol, maxcoeff=maxcoeff, maxsteps=maxsteps)
    if rel is None:
        return None, None, tol
    res = sum(r * v for r, v in zip(rel, vec))
    return rel, res, tol

def show(name, basis, rel, res, tol, D):
    if rel is None:
        print(f"  [{name}] D={D}: NO relation (tol {mp.nstr(tol,2)}, dps {mp.mp.dps})"); return
    mx = max(abs(r) for r in rel)
    print(f"  [{name}] D={D}: max|coeff| = {mx}, residual {mp.nstr(res, 4)} (dps {mp.mp.dps}); relation (target first):")
    print("     ", rel)
    if rel[0] != 0:
        for (nm, _), r in zip(basis, rel[1:]):
            if r != 0:
                print(f"        {nm}: {mp.mpf(-r) / rel[0]}   = {-r}/{rel[0]}")

def level2_by_weight(C):
    """level-2 (alternating MZV) Q-bases by weight 0..6 (weight 6 = basis13)."""
    lam, z3, z5, pi2, Li4, Li5 = (C[x] for x in ("lam", "z3", "z5", "pi2", "Li4", "Li5"))
    one = mp.mpf(1)
    return {0: [("1", one)],
            1: [("lam", lam)],
            2: [("pi^2", pi2), ("lam^2", lam ** 2)],
            3: [("z3", z3), ("pi^2 lam", pi2 * lam), ("lam^3", lam ** 3)],
            4: [("z3 lam", z3 * lam), ("pi^4", pi2 ** 2), ("pi^2 lam^2", pi2 * lam ** 2), ("lam^4", lam ** 4), ("Li4", Li4)],
            5: [("z5", z5), ("pi^2 z3", pi2 * z3), ("z3 lam^2", z3 * lam ** 2), ("pi^4 lam", pi2 ** 2 * lam),
                ("pi^2 lam^3", pi2 * lam ** 3), ("lam^5", lam ** 5), ("lam Li4", lam * Li4), ("Li5", Li5)],
            6: basis13(C)}

def Lfamily(C, Lval, sqrt2=False):
    """L^m * (level-2 weight 6-m), m=0..6, optionally also times sqrt2."""
    B = level2_by_weight(C)
    out = []
    for m in range(0, 7):
        for nm, v in B[6 - m]:
            out.append((f"L^{m}*{nm}", Lval ** m * v))
    if sqrt2:
        r2 = mp.sqrt(2)
        out += [("r2*" + nm, r2 * v) for nm, v in out]
    return out
