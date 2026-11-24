import sys, random
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
from fractions import Fraction as Fr
from cor5_sym import *
import cor4_sym as c4

random.seed(12345)
# --- control 1: 100 random shuffle controls give 0 in delta_22 ---
e1 = (Fr(1), Fr(0)); e2 = (Fr(0), Fr(1))
letters = [e1, e2]
def rw(n):
    return tuple(random.choice(letters) for _ in range(n))
ok = 0
for t in range(100):
    nu = random.choice([1, 2])
    u, v = rw(nu), rw(4 - nu)
    S = sshuffle({u: Fr(random.randint(-3, 3) or 1)}, {v: Fr(random.randint(-3, 3) or 1)})
    S = expand_sym(S, 2)
    d = c4.delta22L2_spec(S, 2)
    if not d:
        ok += 1
print(f"shuffle controls giving delta22 = 0: {ok}/100")

# --- control 2: Li4 at S-unit args gives delta22 = 0 ---
ok2 = 0
import json
sols = json.load(open("sunit23_solutions.json"))
for xs, fx, fy in sols:
    rv = (Fr(fx[1]), Fr(fx[2])); om = (Fr(fy[1]), Fr(fy[2]))
    S = expand_sym(li_sym(4, rv, om), 2)
    if not c4.delta22L2_spec(S, 2):
        ok2 += 1
print(f"Li4(x) at 21 S-unit args giving delta22 = 0: {ok2}/21")
# positive control: delta22 of H_{0,1,0,-1}(-x)-like symbol is nonzero
Z = (Fr(0),)*4
def vec(j, val=1):
    v = list(Z); v[j] = Fr(val); return tuple(v)
G_mc = {0: vec(1), 1: lneg(vec(3)), -1: vec(2)}   # argument -c over basis (2,c,1-c,1+c)
Spos = hpl_sym((0,1,0,-1), G_mc)
d = c4.delta22L2_spec(Spos, 4)
print("positive control delta22(H_{0,1,0,-1}(-c)) nonzero:", bool(d))
