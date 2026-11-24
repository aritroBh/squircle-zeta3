cor4_fix = '''
# ---------- specialisation ----------
def make_letter_map(images):
    """images: list of new-basis vectors = images of the old basis letters.
    Returns function mapping an old letter vector to a new letter vector."""
    def f(u):
        dim = len(images[0])
        return tuple(sum(u[j] * images[j][k] for j in range(len(u))) for k in range(dim))
    return f

def specialize_func(FS, lmap, cpow):
    """FS: function symbol dict (e, word) -> Fr.
    lmap: letter vector -> new letter vector.
    cpow: dict e -> exact field element for c0^e (Fr or (Fr,Fr) pair).
    Returns dict word -> field-element coefficient (same type as cpow values)."""
    out = {}
    for (e, w), cf in FS.items():
        nw = tuple(lmap(v) for v in w)
        c0e = cpow[e]
        if isinstance(c0e, tuple):
            val = (cf * c0e[0], cf * c0e[1])
        else:
            val = cf * c0e
        out[nw] = out.get(nw, 0) + val if not isinstance(val, tuple) else (
            (out.get(nw, (0, 0))[0] + val[0], out.get(nw, (0, 0))[1] + val[1]))
    if out and isinstance(next(iter(out.values())), tuple):
        return {w: c for w, c in out.items() if c != (Fr(0), Fr(0))}
    return {w: c for w, c in out.items() if c != 0}

def delta22_spec(S, dim):
    """(2,2) cobracket for specialised symbols whose coefficients may be (Fr,Fr) pairs."""
    out = {}
    for w, cf in S.items():
        a, b, c, d = w
        for p1, c1 in wedge(a, b, dim).items():
            for p2, c2 in wedge(c, d, dim).items():
                key = (p1, p2)
                if isinstance(cf, tuple):
                    cur = out.get(key, (Fr(0), Fr(0)))
                    out[key] = (cur[0] + cf[0] * c1 * c2, cur[1] + cf[1] * c1 * c2)
                else:
                    out[key] = out.get(key, 0) + cf * c1 * c2
    if out and isinstance(next(iter(out.values())), tuple):
        return {k: v for k, v in out.items() if v != (Fr(0), Fr(0))}
    return {k: v for k, v in out.items() if v}
'''
cor = "<project>/research/08_theory/COR4_M3alpha_classical/"
src = open(cor + "cor4_sym.py").read()
# remove the two broken functions
src = src.split("# ---------- specialisation ----------")[0]
src += cor4_fix
open(cor + "cor4_sym.py", "w").write(src)
import importlib, sys
sys.path.insert(0, cor)
import cor4_sym
importlib.reload(cor4_sym)
print("module OK")
