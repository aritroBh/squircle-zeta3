"""w6 identification step 14: "Gamma-world" tests for c_3, c_2, c_1, c_0 (the leading coefficient c_4 lies in the Q-span of monomials in
{log2, pi^2, zeta3, zeta5}, no Li_n(1/2), no S6).  Basis G_w = monomials of weight w in {lam, pi^2 (wt 2), z3 (wt 3), z5 (wt 5)} and the
families L^m * G_{6-m}  (L = log(1+sqrt2))  with and without sqrt2 multiples (25 / 50 elements).  LLL (flint) at two precisions:
run 400 / D=380 and run 250 / D=235.  Accept only a relation with small coefficients, residual at working precision, identical at both precisions."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
from w6lll import find_rel
def Gworld(lam, pi2, z3, z5):
    one = mp.mpf(1)
    return {0: [("1", one)], 1: [("lam", lam)], 2: [("lam^2", lam ** 2), ("pi^2", pi2)],
            3: [("lam^3", lam ** 3), ("lam pi^2", lam * pi2), ("z3", z3)],
            4: [("lam^4", lam ** 4), ("lam^2 pi^2", lam ** 2 * pi2), ("pi^4", pi2 ** 2), ("lam z3", lam * z3)],
            5: [("lam^5", lam ** 5), ("lam^3 pi^2", lam ** 3 * pi2), ("lam pi^4", lam * pi2 ** 2), ("lam^2 z3", lam ** 2 * z3),
                ("pi^2 z3", pi2 * z3), ("z5", z5)],
            6: [("lam^6", lam ** 6), ("lam^4 pi^2", lam ** 4 * pi2), ("lam^2 pi^4", lam ** 2 * pi2 ** 2), ("pi^6", pi2 ** 3),
                ("lam^3 z3", lam ** 3 * z3), ("lam pi^2 z3", lam * pi2 * z3), ("z3^2", z3 ** 2), ("lam z5", lam * z5)]}
def family(Lval, G, sqrt2):
    out = []
    for m in range(7):
        for nm, v in G[6 - m]:
            out.append((f"L^{m}*{nm}", Lval ** m * v))
    if sqrt2:
        r2 = mp.sqrt(2)
        out += [("r2*" + nm, r2 * v) for nm, v in out]
    return out
for j in (3, 2, 1, 0):
    for sq in (False, True):
        found = []
        for run, D in ((400, 380), (250, 235)):
            mp.mp.dps = D + 20
            lam = mp.log(2); pi2 = mp.pi ** 2; z3 = mp.zeta(3); z5 = mp.zeta(5)
            B = family(mp.log(1 + mp.sqrt(2)), Gworld(lam, pi2, z3, z5), sq)
            cj = mp.mpf(coef(run, j))
            mp.mp.dps = D + 20
            R = find_rel(cj, B, D)
            n = len(B) + 1
            msg = f"c_{j}, run {run}, D={D}, {n} numbers (sqrt2 multiples: {sq}), spurious size ~10^{D/n:.1f}: "
            if not R:
                msg += "no relation"
            else:
                r = R[0]
                msg += f"best relation max|coeff|={r['maxcoeff']} residual {mp.nstr(r['res'],3)} target coeff {r['rel'][0]}"
            print(msg, flush=True)
            found.append(R[0] if R else None)
            if R and R[0]["maxcoeff"] < 10 ** 5:
                print("    TERMS:", R[0]["terms"])
