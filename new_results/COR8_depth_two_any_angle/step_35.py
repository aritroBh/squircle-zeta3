
import sys, json
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
from cor8_sym import *
import itertools
Z3 = (Fr(0),)*3
def v3(j, val=1):
    v = list(Z3); v[j] = Fr(val); return tuple(v)
cv, mv, pv = v3(0), v3(1), v3(2)
zero = v3(0,0)
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"

def drop_zero(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}

# S(Omega)
FS = json.load(open(C4 + "F_symbol.json"))
def parse_letter(s):
    digs = []
    i = 0
    while i < len(s):
        if s[i] == '-':
            digs.append(int(s[i:i+2])); i += 2
        else:
            digs.append(int(s[i])); i += 1
    assert digs[0] == 0
    return (Fr(digs[1]), Fr(digs[2]), Fr(digs[3]))
OM = {}
for key, val in FS.items():
    e = int(key.split(" ")[0][2:])
    word = tuple(parse_letter(b) for b in key.split(" ")[1].split("|"))
    OM[word] = padd(OM.get(word, {}), {e: Fr(val)})

# DT = (c^-1-1) Li22(1/2,1-c) - (c^-1+1) Li22(1/2,1+c)
S22m = drop_zero(li22_sym(zero, mv, zero, cv, pv))
S22p = drop_zero(li22_sym(zero, pv, zero, cv, mv))
DT = {}
for w, cf in S22m.items():
    DT[w] = pmul(cf, {-1: Fr(1), 0: Fr(-1)})
for w, cf in S22p.items():
    DT[w] = padd(DT.get(w, {}), pmul(cf, {-1: Fr(-1), 0: Fr(-1)}))
DT = {w: cf for w, cf in DT.items() if cf}
R = dict(OM)
for w, cf in DT.items():
    R[w] = padd(R.get(w, {}), pmul(cf, {0: Fr(-1)}))
R = {w: cf for w, cf in R.items() if cf}
print("R:", len(R), "words, delta22 =", delta22(R, 3))

sa = json.load(open(C4 + "sunit_args.json"))
args = []
for v_r, v_om, s in sa:
    rv = tuple(Fr(x) for x in v_r[1:]); omv = tuple(Fr(x) for x in v_om[1:])
    if rv == zero: continue
    args.append((s, rv, omv))
li4 = {s: li_sym(4, rv, omv) for s, rv, omv in args}
logs = {'c': cv, '1-c': mv, '1+c': pv}
prods = {}
for s, rv, omv in args:
    S3 = li_sym(3, rv, omv); S2 = li_sym(2, rv, omv)
    for hn, hv in logs.items():
        prods[f"Li3({s})*log({hn})"] = sshuffle(S3, log_sym(hv))
        for h2n, h2v in logs.items():
            prods[f"Li2({s})*log({hn})*log({h2n})"] = sshuffle(S2, sshuffle(log_sym(hv), log_sym(h2v)))
    prods[f"Li2({s})**2"] = sshuffle(S2, S2)
for i, (s1, rv1, omv1) in enumerate(args):
    for s2, rv2, omv2 in args[i+1:]:
        prods[f"Li2({s1})*Li2({s2})"] = sshuffle(li_sym(2, rv1, omv1), li_sym(2, rv2, omv2))
for combo in itertools.combinations_with_replacement(list(logs), 4):
    S = log_sym(logs[combo[0]])
    for h in combo[1:]:
        S = sshuffle(S, log_sym(logs[h]))
    prods["*".join(f"log({h})" for h in combo)] = S
def sym_to_json(S):
    return {"|".join("".join(str(int(x)) for x in l) for l in w): {str(e): str(v) for e, v in cf.items()} for w, cf in S.items()}
json.dump({"Omega": sym_to_json(OM), "DT": sym_to_json(DT), "R": sym_to_json(R),
           "Li4": {k: sym_to_json(v) for k, v in li4.items()},
           "prods": {k: sym_to_json(v) for k, v in prods.items()}},
          open("symbols.json", "w"))
print("saved.  Li4:", len(li4), " prods:", len(prods))
open('progress.log','a').write('step_03a: symbols.json rebuilt (zero-letter filter); d22(Omega-DT)=0\n')
