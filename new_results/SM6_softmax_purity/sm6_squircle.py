# SM6 item 4: same reduction for the p-norm ball.  Exact: I_p = 4 sqrt2 [1 + int_0^inf w(v) ((1+e^{-2v/tau})^tau - 1) dv],
# w(v) = sqrt(1+tanh 2v)/cosh 2v = 2 sqrt2 e^{-2v} (1+e^{-4v})^{-3/2}   (v = artanh(tan phi), phi = pi/4 - theta).
import sympy as sp, pickle
sol, free = pickle.load(open("sm6_reduce.pkl", "rb"))
p = sp.Symbol('p'); v = sp.Symbol('v')
D = 6
g = {(a, b): sp.Symbol('g%d%d' % (a, b)) for a in range(D) for b in range(1, D+1-a)}
w = sp.series(sp.sqrt(1+sp.tanh(2*v))/sp.cosh(2*v), v, 0, 8).removeO()
w2 = sp.series(2*sp.sqrt(2)*sp.exp(-2*v)*(1+sp.exp(-4*v))**sp.Rational(-3, 2), v, 0, 8).removeO()
assert sp.simplify(w-w2) == 0
wk = [sp.nsimplify(sp.expand(w).coeff(v, k)) for k in range(8)]
print("w_k:", wk)
# soft-max weight 1 + v/sqrt(1-v^2)
ws = sp.series(1+v/sp.sqrt(1-v**2), v, 0, 8).removeO()
wsk = [sp.expand(ws).coeff(v, k) for k in range(8)]
print("w_k softmax:", wsk)
def cn(wk, pref, n):
    # I/pref - 1 ~ sum_k w_k k! tau^{k+1} sum_b gamma_kb tau^b ;  coefficient of tau^n
    return sp.expand(pref*sum(wk[k]*sp.factorial(k)*g[k, n-1-k] for k in range(0, n-1) if (k, n-1-k) in g))
Z3, Z5, Z7 = [sp.Symbol('z%d' % i) for i in (3, 5, 7)]
for n in range(2, 8):
    if n-1 > D: break
    e = cn(wk, 4*sp.sqrt(2), n)
    r = sp.simplify(sp.expand(e.subs(sol)))
    es = cn(wsk, 8, n); rs = sp.simplify(sp.expand(es.subs(sol)))
    print("n=%d squircle: %s    free left: %s" % (n, sp.factor(r), r.free_symbols & set(free)))
    print("     softmax: %s    free left: %s" % (rs, rs.free_symbols & set(free)))
