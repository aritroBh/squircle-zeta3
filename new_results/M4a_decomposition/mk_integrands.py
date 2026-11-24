"""Exact integrands for M_k = -F^(k)(0)/(2 k!) in the Euler variable v.

F^(k)(0) = int_0^1 d^k/de^k [a + b - sqrt(a^2+b^2)]_{e=0} dt,  a = e^{(e-1) log t}, b = e^{(e-1) log(1-t)}.
Writing l1 = log t, l2 = log(1-t), the integrand is a polynomial in l1, l2 with algebraic coefficients;
under t = (v^2+2v-1)/(4v) each coefficient times dt/dv is rational in v over Q(sqrt2).
Using t <-> 1-t (v <-> 1/v) every monomial l1^i l2^j with i < j is folded onto l1^j l2^i.

Output: for each k, a dict {(i, j): R_ij(v)} with  M_k = int_b^a sum R_ij(v) l1^i l2^j dv,
together with the partial-fraction structure at v = +-i.
"""
import sympy as sp

e, t, v = sp.symbols('e t v', positive=True)
l1, l2, Asym, Bsym, Qs = sp.symbols('l1 l2 A B Q', positive=True)
S2 = sp.sqrt(2)
a_, b_ = 1 + S2, S2 - 1


def t_integrand(k):
    """d^k/de^k Phi at e=0 as polynomial in l1,l2 with coefficients in t (Q = sqrt(t^2+(1-t)^2))."""
    Phi = Asym * sp.exp(e * l1) + Bsym * sp.exp(e * l2) - sp.sqrt(Asym ** 2 * sp.exp(2 * e * l1) + Bsym ** 2 * sp.exp(2 * e * l2))
    d = sp.diff(Phi, e, k).subs(e, 0)
    d = d.subs(sp.sqrt(Asym ** 2 + Bsym ** 2), Qs / (t * (1 - t)))
    d = d.subs({Asym: 1 / t, Bsym: 1 / (1 - t)})
    d = sp.expand(d)
    poly = sp.Poly(d, l1, l2)
    return {m: sp.simplify(c) for m, c in poly.terms()}


def v_integrand(k):
    """{(i,j): R_ij(v)} folded so that i >= j, for M_k = -F^(k)(0)/(2 k!)."""
    terms = t_integrand(k)
    tv = (v ** 2 + 2 * v - 1) / (4 * v)
    Qv = (v ** 2 + 1) / (2 * S2 * v)
    dtdv = sp.diff(tv, v)
    folded = {}
    for (i, j), c in terms.items():
        if i >= j:
            key, cc = (i, j), c
        else:      # t -> 1-t swaps l1,l2 and t -> 1-t in the coefficient (Q invariant)
            key, cc = (j, i), c.subs(t, 1 - t)
        folded[key] = folded.get(key, 0) + cc
    out = {}
    pref = -sp.Rational(1, 2) / sp.factorial(k)
    for key, c in folded.items():
        R = sp.together(pref * c.subs(Qs, Qv).subs(t, tv) * dtdv)
        R = sp.factor(sp.radsimp(sp.cancel(R)), extension=S2)
        out[key] = R
    return out


if __name__ == "__main__":
    import sys
    for k in [int(x) for x in sys.argv[1:]] or [1, 2, 3]:
        print(f"===== k = {k}")
        Rk = v_integrand(k)
        for key in sorted(Rk, reverse=True):
            R = Rk[key]
            print(f"  l1^{key[0]} l2^{key[1]} :  {R}")
            pf = sp.apart(R, v, extension=[S2, sp.I])
            for term in sp.Add.make_args(pf):
                print("        ", sp.simplify(term))
