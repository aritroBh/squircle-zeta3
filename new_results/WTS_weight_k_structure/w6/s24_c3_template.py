"""w6 identification step 17: c_3 (next-to-leading coefficient of W_6(k)) with the template learned at weight 4:
c_1^(4) = Q-combination of { Gamma-world L-family wt 4 , Li3(x)*{lam, L}, x in {b,-b,1/sqrt2} }  (no Li4, no Li2(b), no sqrt2).
Weight-6 analogue tried here:  Gamma-world L-family wt 6 (25)  +  Li5(x)*{lam, L}, x in {b,-b,1/sqrt2,1/2}  (8)  = 33 numbers, no sqrt2;
then variants (with sqrt2 multiples).  LLL (flint) at D=380 (run 400) and D=235 (run 250).  Printed numbers from this script."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
from w6lll import find_rel
exec(open(os.path.join(W, "s21_gamma_world.py")).read().split("for j in (3, 2, 1, 0):")[0].split("from w6lll import find_rel")[1])   # reuse Gworld(), family()
def run(label, builder, js=(3,), runs=((400, 380), (250, 235))):
    for j in js:
        for run_, D in runs:
            mp.mp.dps = D + 25
            B = builder()
            cj = mp.mpf(coef(run_, j)); mp.mp.dps = D + 25
            n = len(B) + 1
            R = find_rel(cj, B, D)
            msg = f"[{label}] c_{j}, run {run_}, D={D}, {n} numbers, spurious size ~10^{D/n:.1f}: "
            if not R: msg += "no relation"
            else:
                r = R[0]; msg += f"best relation max|coeff|={r['maxcoeff']} residual {mp.nstr(r['res'],3)} target coeff {r['rel'][0]}"
            print(msg, flush=True)
            if R and R[0]["maxcoeff"] < 10 ** 6:
                print("    TERMS:", R[0]["terms"], flush=True)
def build_A(sq=False, extra=False):
    def f():
        lam = mp.log(2); pi2 = mp.pi ** 2; z3 = mp.zeta(3); z5 = mp.zeta(5)
        r2 = mp.sqrt(2); b = r2 - 1; L = mp.log(1 + r2)
        B = family(L, Gworld(lam, pi2, z3, z5), False)
        for nm, x in (("b", b), ("-b", -b), ("1/r2", 1 / r2), ("1/2", mp.mpf(1) / 2)):
            li = mp.polylog(5, x)
            B += [(f"Li5({nm})*lam", li * lam), (f"Li5({nm})*L", li * L)]
        if sq:
            B += [("r2*" + n, r2 * v) for n, v in B]
        return B
    return f
run("A: G-world L-fam + Li5*{lam,L}", build_A(False))
