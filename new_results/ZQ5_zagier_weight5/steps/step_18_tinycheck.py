
import sys; sys.path.insert(0,'.')
from fractions import Fraction
import json
import mpmath as mp
import w5lib as W
import pipe
D=8
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
RM = json.load(open(pipe.cand_dir(D)+"/rel_min.json"))
names = RM["names"]; c = [Fraction(x) for x in RM["c"]]
cmap = dict(zip(names, c))
assert cmap["T"] == -1
d = 2
T = {}
for i, cd in enumerate(cands):
    q = cmap.get(f"Li5({pipe.fmt_qd(cd['x'])})", Fraction(0))
    if not q: continue
    s = W.symbol_li5_dense(cd['vx'], cd['v1m'], d)
    for idx, cf in enumerate(s):
        if cf: T[idx] = T.get(idx, Fraction(0)) + q*cf
T = {k:v for k,v in T.items() if v}
print("(b) exact full symbol of RHS:", "ZERO" if not T else f"NONZERO: {T}")
bad = 0
for i, cd in enumerate(cands):
    s = W.symbol_li5_dense(cd['vx'], cd['v1m'], d)
    dic = {}
    for idx, cf in enumerate(s):
        if cf:
            tup = []
            x = idx
            for _ in range(5): tup.append(x % d); x //= d
            dic[tuple(reversed(tup))] = cf
    for mode in ("full","cut"):
        if W.cobracket32(dic, d, mode=mode): bad += 1
print("(c) Li5 symbols with nonzero (3,2)-cobracket:", bad, "(expect 0: each individually kills delta)")
P = 1300
mp.mp.dps = P+30
lam = mp.log(mp.sqrt(2)); lL = mp.log(1+mp.sqrt(2)); pi2 = mp.pi**2
tot = -(-1)*0 + W.Lchi_hurwitz(D, P)/mp.sqrt(2)   # +T (since c_T=-1 => residual = RHS - T)
for nm, cf in cmap.items():
    if nm == "T" or not cf: continue
    if nm == "z5": v = mp.zeta(5)
    elif nm.startswith("Li5("):
        i = [j for j,cd in enumerate(cands) if f"Li5({pipe.fmt_qd(cd['x'])})" == nm][0]
        v = W.li5_series(cands[i]['x'].emb(P+30), P)
    else:
        v = mp.mpf(1)
        for f in nm.split("*"): v *= {"ls2": lam, "lE": lL, "pi2": pi2}[f]
    tot -= cf*v
print("T - RHS at fresh 1300 digits:", mp.nstr(tot, 4))
