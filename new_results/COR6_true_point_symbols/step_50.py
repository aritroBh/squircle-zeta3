ft = r'''"""cor6_ft.py -- tensors with number-field coefficients, handled componentwise over Q."""
from fractions import Fraction as Fr
from cor6_sym import sadd, sscale, expand_sym, W4Quotient, delta22

def ft_add(*Ts):
    out = {}
    for T in Ts:
        for w, c in T.items():
            cur = out.get(w)
            if cur is None: out[w] = c
            else: out[w] = tuple(x + y for x, y in zip(cur, c))
    return {w: c for w, c in out.items() if any(c)}

def ft_scale(F, a, T):
    """multiply an Fr-tensor (word->Fr) by NF element a -> NF tensor."""
    if all(x == 0 for x in a): return {}
    return {w: tuple(x * cf for x in a) for w, cf in T.items()}

def ft_components(T, d):
    """NF tensor -> list of d Fr tensors (coefficient of theta^j)."""
    comps = [{} for _ in range(d)]
    for w, c in T.items():
        for j in range(d):
            if c[j]:
                comps[j][w] = comps[j].get(w, 0) + c[j]
    return [{w: cf for w, cf in comp.items() if cf} for comp in comps]

def ft_coords(Q, T, d):
    """Q4 coords per theta-power component."""
    return [Q.coords(comp) for comp in ft_components(T, d)]

def ft_delta22(T, d, dim):
    """delta22 per theta-power component; returns dict wedgepair -> NF tuple."""
    out = {}
    for j, comp in enumerate(ft_components(T, d)):
        for k, v in delta22(comp, dim).items():
            cur = out.get(k)
            if cur is None: cur = [Fr(0)] * d
            cur = list(cur); cur[j] += v
            out[k] = tuple(cur)
    return {k: v for k, v in out.items() if any(v)}
'''
with open("cor6_ft.py", "w") as f:
    f.write(ft)

# ---- control (c): square corner c = 1/sqrt2 ----
from fractions import Fraction as Fr
import importlib, json
import cor6_nf, cor6_sunit, cor6_sym, cor6_ft
for m in (cor6_nf, cor6_sunit, cor6_sym, cor6_ft): importlib.reload(m)
from cor6_nf import NF
from cor6_sunit import build_basis, exponent_vector
from cor6_sym import *
from cor6_ft import *

base = "<project>/research/"
D = json.load(open(base + "10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))

F = NF([-2, 0, 1])
d = 2
basis, pool = build_basis(F, [2], 2, verbose=False)
sq2 = F.theta()
c0 = (Fr(0), Fr(1, 2))  # 1/sqrt2
one = F.one()

# letter vector map for this field
letters_needed = {}
def vecF(el):
    key = el
    if key not in letters_needed:
        ev = exponent_vector(F, el, basis)
        assert ev is not None, ("factor failed", el)
        letters_needed[key] = ev
    return letters_needed[key]

# Omega(c0) as NF tensor
SOm = {}
for key, val in D:
    e, w = key.split(' ', 1); e = int(e[2:]); w = eval(w[1:])
    if len(w) != 4: continue
    ce = F.pow(c0, e) if e >= 0 else F.pow(c0, e)
    SOm = ft_add(SOm, ft_scale(F, F.scale(Fr(val), ce), sym_H(w, c0, vecF)))

print("number of words in S(Omega(1/sqrt2)):", len(SOm))
Q2 = W4Quotient(2)
d22 = ft_delta22(SOm, d, 2)
print("delta_22(Omega(1/sqrt2)):", d22, "(Lambda^2 B_2 = 0 in dim 2, so must be empty)")
print("class(Omega(1/sqrt2)) per component:", ft_coords(Q2, SOm, d))
