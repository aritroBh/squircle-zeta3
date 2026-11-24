import sympy as sp
u, S = sp.symbols('u S')  # S = sqrt(1+u^2), S^2 = 1+u^2
S2 = sp.sqrt(2)
a = 1 + S2; b = S2 - 1
v = u + S
def red(expr):
    # reduce expr (rational in u,S over Q(sqrt2)) modulo S^2 = 1+u^2 -> form (P + Q S)/R with P,Q,R polys in u
    expr = sp.together(expr)
    num, den = sp.fraction(expr)
    num = sp.Poly(sp.expand(num), S, extension=S2)
    den = sp.Poly(sp.expand(den), S, extension=S2)
    rel = sp.Poly(S**2 - (1+u**2), S, extension=S2)
    _, numr = sp.div(num, rel)
    _, denr = sp.div(num*0 + den, rel)
    # denr = d0 + d1 S ; invert: 1/(d0+d1 S) = (d0 - d1 S)/(d0^2 - d1^2 (1+u^2))
    d0 = sp.expand(denr.as_expr()).coeff(S, 0)
    d1 = sp.expand(denr.as_expr()).coeff(S, 1)
    norm = sp.expand(d0**2 - d1**2*(1+u**2))
    inv = sp.expand(d0 - d1*S)
    res = sp.expand(numr.as_expr() * inv)
    c0 = sp.simplify(sp.expand(res).coeff(S, 0) / norm)
    c1 = sp.simplify(sp.expand(res).coeff(S, 1) / norm)
    return c0, c1   # measure = c0 + c1 * S
dvdu = v/S  # = u/S + 1
combos = {
 'odd': 1/(v+a) - 1/(v-a),
 'even': 1/(v+a) - 1/(v+b),
 'diag': 1/(v+a) - 1/(2*v),
 'oddlast': 1/(v+a) - 1/(2*(v-a)) + 1/(2*(v+b)) - 1/(2*v),
}
for name, Rv in combos.items():
    c0, c1 = red(Rv * dvdu)
    # target form: (A + B u + (C + D u)/S)/(1-u^2) ; i.e. c0*(1-u^2) = A + B u + C u? no:
    # c0 + c1 S = [A + B u + (C + D u)/S]/(1-u^2) = [A + B u]/(1-u^2) + [(C+Du)]/[(1-u^2) S]
    # multiply: c0 = (A+Bu)/(1-u^2) + C' /(1-u^2)?? better: write c0 = (A + B u + E*(1+u^2)? ) ... solve:
    # c0*(1-u^2) = A + B u + C u? no. c0*(1-u^2) should equal A + B u + 0/S ... and c1*(1-u^2) = (C + D u)/(1+u^2)?? 
    # simpler: (c0 + c1 S)(1-u^2) = A + Bu + (C+Du)/S  =>  c0(1-u^2) = A + Bu + C/S*? no.
    # (C+Du)/S = (C+Du) S/(1+u^2). So rational part: A + Bu = c0(1-u^2) - (C+Du)*0 ... do:
    # c0(1-u^2) = A + Bu  (poly in u, deg<=1)
    # c1(1-u^2) = (C + D u)/(1+u^2)
    A_B = sp.expand(c0*(1-u**2))
    print(name, ': c0*(1-u^2) =', A_B, '   c1*(1-u^2)*(1+u^2) =', sp.expand(c1*(1-u**2)*(1+u**2)))
