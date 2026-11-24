"""cor6_ft.py -- tensors with number-field coefficients, handled componentwise over Q."""
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
