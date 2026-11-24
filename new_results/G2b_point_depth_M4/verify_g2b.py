"""verify_g2b.py -- < 5 min re-verification of the G2b results.

1. CONTROL: the fixed point-level (3,2) pairing of Omega_5(1/sqrt(2)) == (84, -38 sqrt(2))
   (M4c's value), and the pairing code on M4c's own symbol.json == same.
2. Pairings at c0 = 1/2 (hexagon) == (41, -33/2) and c0 = 1/3 == (1838/3, -284).
3. Shuffle and Li5 controls pair to zero.
4. N5 cobracket match: N5 = (5/4)*Li32(1/4,4/3) + (1/3)*Li32(8/9,1) has
   quotient cobracket (-41, -33/2) and symbol pairing (41, -33/2) (both exact).
Prints PASS/FAIL per item; exits nonzero on any FAIL.
"""
import sys, os, json, time
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
G2B = BASE + "/G2b_point_depth_M4/"
sys.path.insert(0, BASE + "/COR6_true_point_symbols")
sys.path.insert(0, BASE + "/G2_M4alpha")
os.chdir(BASE + "/G2_M4alpha")
import importlib.util
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); mod.__name__ = name
    spec.loader.exec_module(mod)
    return mod
s01 = load("g2b_step01", G2B + "step_01.py")
s05 = load("g2b_step05", G2B + "step_05.py")
from step_13 import W5Point
from step_12 import br

t0 = time.time()
ok_all = True
def check(name, cond):
    global ok_all
    ok_all = ok_all and bool(cond)
    print(f"  [{'PASS' if cond else 'FAIL'}] {name} ({time.time()-t0:.0f}s)", flush=True)

print("verify_g2b:")

# ---- 1. control at c0 = 1/sqrt(2)
sj = json.load(open(BASE + "/M4c_w5_closed_form/symbol.json"))
Sm4c = {tuple(int(ch) for ch in w): (Fr(p), Fr(q)) for w, (p, q) in sj["symbol"].items()}
o_m4c = s01.obstruction_fixed(Sm4c)
check("pairing code on M4c symbol.json == (84, -38 sqrt2)",
      o_m4c[0] == (Fr(84), Fr(0)) and o_m4c[1] == (Fr(0), Fr(-38)))

P4 = W5Point("c=1/sqrt2", [-2, 0, 1], (Fr(0), Fr(1, 2)),
             basis=[(Fr(2), Fr(0)), (Fr(1), Fr(1))])
obs4 = s01.obstruction_fixed(s01.ft_expand_sym(P4.omega5(), 2))
check("CONTROL: fixed point pairing at 1/sqrt(2) == (84, -38 sqrt2)",
      obs4[0] == (Fr(84), Fr(0)) and obs4[1] == (Fr(0), Fr(-38)))

# ---- 2. hexagon and c = 1/3 pairings
P6 = W5Point("c=1/2", [-1, 1], (Fr(1, 2),))
obs6 = s01.obstruction_fixed(s01.ft_expand_sym(P6.omega5(), 2))
check("hexagon c=1/2 pairing == (41, -33/2)",
      obs6[0] == (Fr(41),) and obs6[1] == (Fr(-33, 2),))
P3 = W5Point("c=1/3", [-1, 1], (Fr(1, 3),))
obs3 = s01.obstruction_fixed(s01.ft_expand_sym(P3.omega5(), 2))
check("c=1/3 pairing == (1838/3, -284)",
      obs3[0] == (Fr(1838, 3),) and obs3[1] == (Fr(-284),))

# ---- 3. shuffle and Li5 controls (reduced counts)
check("shuffle control (60 random products)", s01.shuffle_control(60))
check("Li5 control (60 random Li5-form symbols)", s01.li5_control(60))

# ---- 4. N5 cobracket match (exact)
X = {(0,): Fr(1)}; Y = {(1,): Fr(1)}
P = br(X, Y); Q = br(X, br(X, Y)); T = br(Y, br(X, Y))
def obs_of_SE(SE):
    SE = {w: (cf[0] if isinstance(cf, tuple) else cf) for w, cf in SE.items()}
    out = []
    for A in (Q, T):
        tot = Fr(0)
        for u, cu in A.items():
            for v, dv in P.items():
                tot += cu * dv * (SE.get(u + v, Fr(0)) - SE.get(v + u, Fr(0)))
        out.append(tot)
    return out
N5 = [(Fr(1, 4), Fr(4, 3), Fr(5, 4)), (Fr(8, 9), Fr(1), Fr(1, 3))]
qtot = [Fr(0), Fr(0)]; otot = [Fr(0), Fr(0)]
for x, y, d in N5:
    dw = s05.delta_word_sym(s05.li32_word(x, y), 0, 1)
    qcw = s05.quotient_coords(s05.ft_expand_sym({w: (cf,) for w, cf in dw.items()}, 2))[0]
    S_at = s05.ft_expand_sym({w: (cf,) for w, cf in s05.li32_sym(x, y).items()}, 2)
    ob = obs_of_SE(S_at)
    qtot = [qtot[0] + d * qcw[0], qtot[1] + d * qcw[1]]
    otot = [otot[0] + d * ob[0], otot[1] + d * ob[1]]
check("N5 quotient cobracket == (-41, -33/2) (== Delta(R4^(5)))",
      tuple(qtot) == (Fr(-41), Fr(-33, 2)))
check("N5 symbol pairing == (41, -33/2) (== obs(Omega_5(1/2)))",
      tuple(otot) == (Fr(41), Fr(-33, 2)))
check("N5 atoms convergent in the paper convention (|x|<1, |xy|<1)",
      all(abs(x) < 1 and abs(x * y) < 1 for x, y, d in N5))

print(f"verify_g2b: {'ALL PASS' if ok_all else 'FAILURES'} ({time.time()-t0:.0f}s)")
sys.exit(0 if ok_all else 1)
