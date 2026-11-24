import sympy as sp
u = sp.symbols('u')
S2 = sp.sqrt(2)
a = 1 + S2; b = S2 - 1
def mul(P, Q):
    p,q = P; r,t = Q
    return (sp.expand(p*r + q*t*(1+u**2)), sp.expand(p*t + q*r))
def inv(P):
    p,q = P
    norm = sp.expand(p**2 - q**2*(1+u**2))
    return (p/norm, -q/norm)
def div(P, Q): return mul(P, inv(Q))
def add(P, Q): return (sp.together(P[0]+Q[0]), sp.together(P[1]+Q[1]))
def sub(P, Q): return (sp.together(P[0]-Q[0]), sp.together(P[1]-Q[1]))
def scalar(r): return (sp.sympify(r), sp.Integer(0))
v = (u, sp.Integer(1))  # u + S
dvdu = div(v, (sp.Integer(0), sp.Integer(1)))  # v/S
def pole(r):
    # 1/(v - r) as pair
    return inv(sub(v, scalar(r)))
combos = {
 'odd': sub(pole(-a), pole(a)),
 'even': sub(pole(-a), pole(-b)),
 'diag': sub(pole(-a), div(scalar(sp.Rational(1,2)), v)),
 'oddlast': sub(add(sub(pole(-a), div(scalar(sp.Rational(1,2)), sub(v, scalar(a)))),
                  div(scalar(sp.Rational(-1,2)), v)),
               div(scalar(sp.Rational(-1,2)), add(v, scalar(b)))),
}
for name, R in combos.items():
    M = mul(R, dvdu)
    p, q = sp.simplify(M[0]), sp.simplify(M[1])
    # check form: p*(1-u^2) poly deg<=1, q*(1-u^2)*(1+u^2) poly deg<=1
    P = sp.expand(p*(1-u**2)); Q = sp.expand(q*(1-u**2)*(1+u**2))
    print(name, ': [A+Bu] =', sp.expand(P), '   [C+Du] =', sp.expand(Q))
