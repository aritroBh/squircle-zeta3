# Planner check of E5: P^(k)(1) for K_p = A·{|x|^p+|y|^p<=1}, straight from the arc-length integral and log-Beta moments.
import sympy as sp
s, X, Y, e1, e2, d = sp.symbols('s X Y e1 e2 d')
z2, z3, z5 = sp.symbols('z2 z3 z5')
z4 = sp.Rational(2,5)*z2**2
N = 5
zeta = {2: z2, 3: z3, 4: z4, 5: z5}
lg = lambda e: sum((-1)**k*zeta[k]*e**k/k for k in range(2, N+1))   # log Gamma(1+e) + gamma e
B = sp.exp(lg(e1)+lg(e2)-lg(e1+e2))/(1+e1+e2)
t = sp.symbols('t')
Bs = sp.expand(sp.series(B.subs({e1: t*e1, e2: t*e2}), t, 0, N+1).removeO().subs(t, 1))
def M(i, j):  # int_0^1 log^i t log^j (1-t) dt
    return sp.factorial(i)*sp.factorial(j)*Bs.coeff(e1, i).coeff(e2, j) if i+j <= N else None
def perimeter_series(A):
    G = A.T*A
    total = 0
    for s1 in (1, -1):
        for s2 in (1, -1):
            sig = s1*s2
            F = G[0,0]*sp.exp(2*s*X) - 2*sig*G[0,1]*sp.exp(s*(X+Y)) + G[1,1]*sp.exp(2*s*Y)
            ser = sp.expand(sp.series(sp.sqrt(F), s, 0, N+1).removeO())
            # integrate each monomial in X, Y
            poly = sp.Poly(ser, X, Y, s)
            val = 0
            for (i, j, k), c in poly.terms():
                val += c*M(i, j)*s**k
            total += val
    # P(p) = (1/p) * total, s = 1/p - 1, p = 1 + d
    P = sp.expand(sp.series((total/(1+d)).subs(s, 1/(1+d)-1), d, 0, N+1).removeO())
    return P
for name, A in [("diamond", sp.eye(2)), ("sheared (2,0),(1,2)", sp.Matrix([[2, 1], [0, 2]]))]:
    P = perimeter_series(A)
    print("==", name)
    for k in (3, 5):
        Dk = sp.expand(sp.factorial(k)*P.coeff(d, k))
        Dk = sp.expand(Dk)
        print(f" P^({k})(1):")
        for mono in ([z5, z2*z3, z2**2, z3, z2, 1] if k == 5 else [z3, z2, 1]):
            if mono == 1:
                c = Dk.subs({z2: 0, z3: 0, z5: 0})
            else:
                c = sp.Poly(Dk, z2, z3, z5).coeff_monomial(mono)
            print("   ", mono, ":", sp.nsimplify(sp.radsimp(sp.simplify(c))))
