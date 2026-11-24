"""Depth-two audit, part 2: regenerate the four atoms' pairing vectors from scratch with the W4-session
symbol engine (A3_W4_proof/w4_symbol.py, read-only import) and MY OWN pairing code (not WT6S's, not M5a's w6_symbol.py).
Convention (M5a): Li_{m,n}(x,y) = G_0(0^{m-1}, 1/x, 0^{n-1}, 1/(xy); 1).
Pairings: (3,3) = <S, [X,P] (x) [Y,P] - [Y,P] (x) [X,P]>; (4,2) = <S, C (x) P - P (x) C>, C in {[X,[X,P]], [Y,[X,P]], [Y,[Y,P]]}.
"""
import sys, json, random
from fractions import Fraction as F
sys.dont_write_bytecode = True
TH = "<project>/research/08_theory/"
sys.path.insert(0, TH + "A3_W4_proof")
from w4_symbol import Q, R, A, B, ZERO, ONE, add, scale, shuffle, log_symbol, concat, g_symbol, li_symbol

X = {(0,): ONE}; Y = {(1,): ONE}
def br(a, b): return add(concat(a, b), scale(concat(b, a), -1))
P = br(X, Y)
XP = br(X, P); YP = br(Y, P)
C42 = [br(X, XP), br(Y, XP), br(Y, YP)]   # a, b, c  ([Y,[X,P]] as in M5a)

def pair(S, U, V):
    tot = ZERO
    for u, a in U.items():
        for v, b in V.items():
            tot = tot + a * b * (S.get(u + v, ZERO) - S.get(v + u, ZERO))
    return tot

def pv(S):
    return [pair(S, XP, YP)] + [pair(S, C, P) for C in C42]

def fmt(v): return "%s" % (v if v.q == 0 else (f"({v.q})*sqrt2" if v.p == 0 else f"({v.p}+{v.q}*sqrt2)")) if hasattr(v, "p") else str(v)
def ev(v): return "[" + ", ".join(fmt(x) for x in v) + "]"

def li_mn(m, n, x, y):
    w = (ZERO,) * (m - 1) + (ONE / x,) + (ZERO,) * (n - 1) + (ONE / (x * y),)
    return g_symbol(ZERO, w, ONE)

# controls
random.seed(20261028)
def rnd(n): return {tuple(random.randint(0, 1) for _ in range(n)): Q(F(random.randint(-5, 5))) for _ in range(4)}
bad = 0
for _ in range(300):
    a = random.randint(1, 5); s = shuffle(rnd(a), rnd(6 - a)); bad += any(v != ZERO for v in pv(s))
print("control: shuffle products with nonzero pairing:", bad, "/ 300  (expect 0)")
units = [R - 1, 1 - R, 3 - 2 * R, R / 2, Q(F(1, 2)), 1 - R / 2, 2 - R, Q(F(1, 2)) - R / 4, -2 + 2 * R, 3 + 2 * R]
print("control: Li6(S-unit) with nonzero pairing:", sum(any(v != ZERO for v in pv(li_symbol(6, z))) for z in units), "/", len(units), " (expect 0)")
b = R - 1
print("control: Li42(b,1) =", ev(pv(li_mn(4, 2, b, ONE))), " expect [0,0,0,-1/4]")
print("control: Li33(b,1) =", ev(pv(li_mn(3, 3, b, ONE))), " expect [0,0,0,3/4]")

atoms = {
    "Li33(b,1)": (3, 3, b, ONE),
    "Li33(1/sqrt2,1)": (3, 3, R / 2, ONE),
    "Li33(b,2)": (3, 3, b, Q(2)),
    "Li42(b,1/sqrt2)": (4, 2, b, R / 2),
}
m5a = {
    "Li33(b,1)": ["0", "0", "0", "3/4"],
    "Li33(1/sqrt2,1)": ["0", "3/16", "0", "0"],
    "Li33(b,2)": ["1", "0", "0", "-3"],
    "Li42(b,1/sqrt2)": ["1/4", "-1/8", "1/8", "3/4"],
}
allok = True
out = {}
for name, (m, n, x, y) in atoms.items():
    v = pv(li_mn(m, n, x, y))
    want = [Q(F(s)) for s in m5a[name]]
    ok = all(v[i] == want[i] for i in range(4))
    allok &= ok
    out[name] = [str(fmt(t)) for t in v]
    print(f"{name:18s} regenerated {ev(v)}  M5a table {m5a[name]}  match={ok}")
print("ALL FOUR ATOM PAIRING VECTORS REGENERATED INDEPENDENTLY AND MATCH M5a TABLE:", allok)
json.dump(out, open("audit2_atoms_regenerated.json", "w"), indent=1)
