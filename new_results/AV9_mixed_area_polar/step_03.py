
# step_03.py : EXACT self-dual eta-expansion coefficients c_n (fast manual series arithmetic)
import sympy as sp
import pickle, time
t0=time.time()
N = 8
eta, t, x, y = sp.symbols('eta t x y')
lv, lw = sp.symbols('lv lw')
L2 = sp.log(2)
zeta = sp.zeta

def tscale(poly, deg):
    return sp.expand(poly.subs({x:t*x, y:t*y}))

def series_exp_t(poly_in_t, deg):
    # poly_in_t: polynomial in t (coeffs may be polynomials in other vars), zero constant term
    # returns exp(poly_in_t) truncated to t^deg
    e = [sp.Integer(1)]
    pc = [sp.expand(poly_in_t).coeff(t, k) for k in range(1, deg+1)]
    for n in range(1, deg+1):
        en = sp.expand(sum(k*pc[k-1]*e[n-k] for k in range(1, n+1)))/n
        e.append(sp.expand(en))
    return sum(e[k]*t**k for k in range(deg+1))

def series_log1p_t(poly_in_t, deg):
    # log(1+u), u has zero constant term in t
    res = sp.Integer(0)
    um = poly_in_t
    for m in range(1, deg+1):
        res += sp.Rational((-1)**(m+1), m)*um
        res = sp.expand(res)
        # truncate
        res = sum(sp.expand(res).coeff(t,k)*t**k for k in range(deg+1))
        um = sp.expand(um*poly_in_t)
        um = sum(um.coeff(t,k)*t**k for k in range(deg+1))
    return res

# ---- Beta(1,1) integrals
maxdeg = N
logkern = sum((-1)**k * zeta(k)/k * (x**k + y**k - (x+y)**k) for k in range(2, maxdeg+1))
logkern_t = tscale(logkern, maxdeg)
exppart = series_exp_t(logkern_t, maxdeg)
geom = sum((-1)**m * t**m*(x+y)**m for m in range(maxdeg+1))
Bfull = sp.expand(exppart*geom)
Bfull = sum(Bfull.coeff(t,k)*t**k for k in range(maxdeg+1)).subs(t,1)
Bpoly = sp.Poly(sp.expand(Bfull), x, y)
int_tab = {}
for (a,b), coeff in Bpoly.terms():
    int_tab[(a,b)] = coeff*sp.factorial(a)*sp.factorial(b)
print("beta integrals done", time.time()-t0, flush=True)

# ---- G(eta) = log(e^{eta lv}+e^{eta lw}) = log2 + log(1+u), u = sum p_k eta^k/(2 k!)
u_eta = sum((lv**k+lw**k)/sp.factorial(2*k) * eta**k for k in range(1, N+1))  # WRONG factor, fix below
u_eta = sum((lv**k+lw**k)*sp.Rational(1,2)/sp.factorial(k) * eta**k for k in range(1, N+1))
# log(1+u) via substitution u -> t u then series_log1p in t
u_t = sp.expand(u_eta.subs({lv: lv*t, lw: lw*t}))  # WRONG scaling; must scale eta
# correct: introduce t multiplying eta
u_t = sp.expand(u_eta.subs(eta, eta*t))
Gtail = series_log1p_t(u_t, N).subs(t, 1)
Gser = sp.expand(L2 + Gtail)
g = {0: L2}
for j in range(1, N+1):
    g[j] = sp.expand(Gser).coeff(eta, j)
print("g_j done", time.time()-t0, flush=True)

A = lv + lw
h = {}
for n in range(1, N+1):
    hn = sum((-1)**(n-j)*g[j] for j in range(0, n+1)) + (-1)**n*A + sp.Rational((-1)**n, n)
    h[n] = sp.expand(hn)
print("h_n done", time.time()-t0, flush=True)

# ---- exp of Hser in eta via recurrence
Hcoeffs = [h[n] for n in range(1, N+1)]
e = [sp.Integer(1)]
for n in range(1, N+1):
    en = sp.expand(sum(k*Hcoeffs[k-1]*e[n-k] for k in range(1, n+1)))/n
    e.append(sp.expand(en))
print("exp done", time.time()-t0, flush=True)

c_exact = {}
for n in range(1, N+1):
    poly = sp.Poly(sp.expand(e[n]), lv, lw)
    val = 0
    for (a,b), coeff in poly.terms():
        val += coeff*int_tab.get((a,b), 0)
    c_exact[n] = sp.expand(4*val)
    print(f"c_{n} done", time.time()-t0, flush=True)
    print(f"c_{n} =", c_exact[n], flush=True)
    print("   numeric =", sp.N(c_exact[n], 30), flush=True)

pickle.dump({n: str(c_exact[n]) for n in c_exact}, open('scratchpad/c_exact_selfdual.pkl','wb'))
print("TOTAL", time.time()-t0)
