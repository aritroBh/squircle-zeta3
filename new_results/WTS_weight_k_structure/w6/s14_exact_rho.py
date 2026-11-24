"""w6 identification step 7: exact polynomials in k of every kernel residue of K_50, K_41, K_32 (k >= 6), by exact rational
interpolation on k=6..10 and exact verification on k=11..15.  W_6(k) = sum rho_{key,letter,comp}(k) * N_{key,letter,comp}
with N_{key,letter,0} = int_b^a dv/(v-r) log^i t log^j(1-t) (pole letter r) and the p1 component multiplying sqrt2*N.
Writes rho_exact.json (strings of Fractions, low-to-high degree in k).  No floating point."""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
KEYS = ["5,0", "4,1", "3,2"]
LET = {"0": (0, 0), "a": (1, 1), "-a": (-1, -1), "b": (-1, 1), "-b": (1, -1)}
inv = {v: k for k, v in LET.items()}
def table(k):
    d = json.load(open(os.path.join(W, f"finals_{k}.json")))["finals"]
    out = {}
    for key in KEYS:
        P = d[key]
        assert not P["ipos"] and not P["ineg"] and all(Fr(x) == 0 for x in P["const"])
        for r, c in P["real"]:
            assert Fr(r[2]) == 0 and Fr(r[3]) == 0 and Fr(c[2]) == 0 and Fr(c[3]) == 0
            out[(key, inv[(int(Fr(r[0])), int(Fr(r[1])))])] = (Fr(c[0]), Fr(c[1]))
    return out
T = {k: table(k) for k in range(6, 16)}
def interp(pts):
    n = len(pts); co = [Fr(0)] * n
    for a in range(n):
        ka, ya = pts[a]; num = [Fr(1)]; den = Fr(1)
        for c in range(n):
            if c == a: continue
            kc = pts[c][0]; new = [Fr(0)] * (len(num) + 1)
            for t, v in enumerate(num): new[t + 1] += v; new[t] -= kc * v
            num = new; den *= (ka - kc)
        for t, v in enumerate(num): co[t] += ya * v / den
    return co
rho = {}
for key in KEYS:
    for let in LET:
        for comp in (0, 1):
            ys = {k: T[k].get((key, let), (Fr(0), Fr(0)))[comp] for k in range(6, 16)}
            if all(v == 0 for v in ys.values()): continue
            c = interp([(Fr(k), ys[k]) for k in range(6, 11)])
            for k in range(11, 16):
                assert sum(c[j] * Fr(k) ** j for j in range(5)) == ys[k]      # exact check
            rho[f"{key}|{let}|{comp}"] = [str(x) for x in c]
            print(f"({key}) pole {let:>2} comp p{comp}: " + "  ".join(f"k^{j}:{c[j]}" for j in range(5)))
json.dump(rho, open(os.path.join(W, "rho_exact.json"), "w"), indent=1)
print("all 15 residue polynomials verified exactly on k=11..15")
