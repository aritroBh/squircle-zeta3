# Referee items 2 + 3: numerical spot checks of the window/outside bounds (own code).
from mpmath import mp, mpf, cos, sin, tan, atanh, sec, log, log1p, exp, pi, sqrt, nstr, norm
mp.dps = 120
# ---- squircle
print("== squircle ==")
# h_p = ||n||_q vs direct sup over the p-ball boundary (coarse check)
for p in (mpf('1.3'), mpf('2.5')):
    q = p/(p-1); th = mpf('0.4'); n = (cos(th), sin(th))
    hq = (abs(n[0])**q+abs(n[1])**q)**(1/q)
    best = max((cos(a)/ (abs(cos(a))**p+abs(sin(a))**p)**(1/p))*n[0] + (sin(a)/(abs(cos(a))**p+abs(sin(a))**p)**(1/p))*n[1]
               for a in [2*pi*k/20000 for k in range(20000)])
    print(f" p={p}: ||n||_q={nstr(hq,15)}  grid sup over boundary={nstr(best,15)}")
K = (sqrt(2)-1-pi**2/32)/(pi/4)**4
print(" K =", nstr(K,10))
worst = 0; worstK = 0
for i in range(1,401):
    x = pi/8*i/400
    Lam = atanh(tan(x))
    # Lambda' = sec 2x check by finite difference
    worst = max(worst, (Lam-x)/((2*sqrt(2)/3)*x**3))
    worstK = max(worstK, (Lam - x - 2*x**3/3)/((16*K/5)*x**5))
    assert Lam - x >= 0
print(" max (Lam-x)/((2sqrt2/3)x^3) on (0,pi/8] =", nstr(worst,6), "(must be <=1)")
print(" max (Lam-x-2x^3/3)/((16K/5)x^5)         =", nstr(worstK,6), "(must be <=1)")
print(" log r + 2 atanh(tan x) at x=0.3:", nstr(log(tan(pi/4-mpf('0.3')))+2*atanh(tan(mpf('0.3'))),3))
print(" Lambda'(0.3) - sec(0.6):", nstr(mp.diff(lambda x: atanh(tan(x)), mpf('0.3'))-sec(mpf('0.6')),3))
def L_sq(s,tau):
    x = tau*abs(s); r = tan(pi/4-x); return tau*log1p(r**(1/tau))
g1 = lambda s: log1p(exp(-2*abs(s)))
g3sq = lambda s: -mpf(4)/3*abs(s)**3/(1+exp(2*abs(s)))
mx3 = 0; mx4 = 0; mx3b=0
for tau in (mpf('0.1'), mpf('0.01'), mpf('0.001'), mpf('1e-4')):
    for j in range(1,300):
        s = mpf(j)/10
        if tau*s > pi/8: break
        L = L_sq(s,tau)
        e3 = abs(L - tau*g1(s))
        mx3b = max(mx3b, e3/((4*sqrt(2)/3)*tau**3*s**3*exp(-2*s)))
        mx3 = max(mx3, e3/(tau**3*exp(-s)))
        mx4 = max(mx4, abs(L - tau*g1(s) - tau**3*g3sq(s))/(tau**5*(1+s)**6*exp(-2*s)))
print(" N=3: max |L-tau g1| / ((4sqrt2/3) tau^3 s^3 e^{-2s}) =", nstr(mx3b,6), "(<=1 claimed)")
print(" N=3: max |L-tau g1| / (tau^3 e^{-s})                  =", nstr(mx3,6), f"(claimed constant (4sqrt2/3)(27/e^3)={nstr(4*sqrt(2)/3*27/exp(3),6)})")
print(" N=4: max |L-tau g1-tau^3 g3| / (tau^5 (1+s)^6 e^{-2s}) =", nstr(mx4,6), "(bounded => O(tau^5) claim consistent)")
# outside windows: |theta|<=pi/8, 0<=h-m<=tau*(sqrt2-1)^{1/tau}
mo = 0
for tau in (mpf('0.5'), mpf('0.1'), mpf('0.02')):
    for i in range(0,201):
        th = pi/8*i/200; q = 1/tau
        hp = (cos(th)**q+sin(th)**q)**tau; m = cos(th)
        assert hp - m >= 0
        mo = max(mo, (hp-m)/(tau*(sqrt(2)-1)**q))
print(" outside: max (h-m)/(tau (sqrt2-1)^{1/tau}) =", nstr(mo,6), "(<=1 claimed)")
# ---- soft-max
print("== soft-max ==")
V = [(1,1),(1,-1),(-1,1),(-1,-1)]
def h_sm(th,tau):
    d = [v[0]*cos(th)+v[1]*sin(th) for v in V]; m = max(d)
    return m*sum(exp((x-m)/tau) for x in d)**tau, m
def L_sm(s,tau):
    hh,m = h_sm(tau*s,tau); return log(hh/m)
def L_formula(s,tau):
    u = tau*s; b = exp(-2*cos(u)/tau)*(1+exp(-2*abs(sin(u))/tau))
    return tau*log(1+exp(-2*sin(tau*abs(s))/tau)+b)
md = 0
for tau in (mpf('0.3'), mpf('0.05')):
    for j in range(-50,51):
        s = mpf(j)/mpf(50)*pi/(4*tau)
        md = max(md, abs(L_sm(s,tau)-L_formula(s,tau)))
print(" window formula: max |L(def) - L(formula)| =", nstr(md,3))
g3sm = lambda s: abs(s)**3/(3*(1+exp(2*abs(s))))
m3=0; m4=0
for tau in (mpf('0.1'), mpf('0.01'), mpf('0.001')):
    for j in range(1,400):
        s = mpf(j)/10
        if tau*s > pi/4: break
        L = L_formula(s,tau)
        m3 = max(m3, abs(L-tau*g1(s))/((tau**3*s**3/3)*exp(-mpf('1.79')*s) + 2*tau*exp(-sqrt(2)/tau)))
        m4 = max(m4, abs(L-tau*g1(s)-tau**3*g3sm(s))/(tau**5*(1+s)**6*exp(-mpf('1.79')*s) + tau*exp(-sqrt(2)/tau)))
print(" N=3: max |L-tau g1| / (claimed bound) =", nstr(m3,6), "(<=1 claimed)")
print(" N=4: max |L-tau g1-tau^3 g3| / (tau^5(1+s)^6 e^{-1.79s} + tau e^{-sqrt2/tau}) =", nstr(m4,6), "(bounded)")
print(" (pi/4)^2/3 =", nstr((pi/4)**2/3,6), "  2-(pi/4)^2/3 =", nstr(2-(pi/4)**2/3,6))
