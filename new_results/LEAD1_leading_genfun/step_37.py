import sympy as sp
u = sp.symbols('u')
s = sp.sqrt(1+u**2)
v = u + s
S2 = sp.sqrt(2)
a = 1 + S2; b = S2 - 1
dvdu = sp.simplify(v/s)   # dv/du
print("dv/du =", dvdu)
def measure(Rv):
    return sp.simplify(sp.expand(Rv * dvdu))
combos = {
 'odd (j<(d-1)/2): 1/(v+a)-1/(v-a)': 1/(v+a) - 1/(v-a),
 'even (j<d/2): 1/(v+a)-1/(v+b)': 1/(v+a) - 1/(v+b),
 'diag (even): 1/(v+a)-1/(2v)': 1/(v+a) - 1/(2*v),
 'oddlast: 1/(v+a)-1/(2(v-a))+1/(2(v+b))-1/(2v)': 1/(v+a) - 1/(2*(v-a)) + 1/(2*(v+b)) - 1/(2*v),
}
for name, Rv in combos.items():
    m = sp.simplify(sp.expand(measure(Rv), s))
    # normalize s = sqrt(1+u^2); express as (A + B u + (C + D u)/s)/(1-u^2) with A,B,C,D in Q(sqrt2)
    m2 = sp.radsimp(sp.together(m))
    print(name)
    print('   measure du-factor =', m2)
