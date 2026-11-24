# R100: mixed area V(K_p, K_p') of the squircle and its polar.
# Own derivation: V(K,L) = (1/2) int h_L dS_K, h_{K_p'}(u) = ||u||_p; on the arc x=t^{1/p}, y=(1-t)^{1/p}
# the normal vector n dt has ||n||_p = (1/p)[t^{1-p}+(1-t)^{1-p}]^{1/p}, so
#   V(p) = (2/p) int_0^1 [t^{1-p} + (1-t)^{1-p}]^{1/p} dt.
from mpmath import mp, mpf, quad, pi, log, zeta, exp, fsum, nstr, mpc
mp.dps = 50
# By symmetry t <-> 1-t and t = exp(-x):  V(p) = (4/p) int_{log 2}^inf e^{-x/p} (1 + (e^x - 1)^{1-p})^{1/p} dx,
# (cut at x = 130|p|: the tail is below e^{-125}) a smooth, exponentially decaying integrand (the raw t-form loses digits near t = 0 when p is large).
def V(p):
    l2 = log(2)
    X = 130*max(1, abs(p))  # decay rate is Re(1/p)
    return 4/p*quad(lambda x: exp(-x/p)*(1+(exp(x)-1)**(1-p))**(1/p), [l2, l2+1, l2+4, l2+16, l2+48, l2+128, X])
def V_eta(eta):  # V is symmetric in p, p'; use q = p' = 1+eta.
    # For complex eta the integrand has branch points where (e^x-1)^(-eta) = -1, i.e. |eta| ~ pi/x; on |eta| = 0.02
    # they sit beyond x ~ 157, so we cut at x = 130 (tail < e^{-127}).
    q = 1+eta; l2 = log(2)
    return 4/q*quad(lambda x: exp(-x/q)*(1+(exp(x)-1)**(1-q))**(1/q), [l2, l2+1, l2+4, l2+16, l2+48, 130])
print("V(2) - pi =", nstr(V(mpf(2))-pi, 3), "  V(p=1) =", nstr(V(mpf(1)), 20))
for eta in [mpf('0.1'), mpf('0.3'), mpf('0.5')]:
    p = 1+1/eta
    print("duality eta=%s: V(p)=%s  V(p')=%s  diff %s" % (eta, nstr(V(p), 25), nstr(V(1+eta), 25), nstr(V(p)-V(1+eta), 3)))
def cauchy(f, r, N, K):
    vals = [f(r*exp(2j*pi*j/N)) for j in range(N//2+1)]
    full = vals + [mpc(v).conjugate() for v in vals[1:N//2][::-1]]
    return [(fsum(full[j]*exp(-2j*pi*j*k/N) for j in range(N))/N/r**k).real for k in range(K+1)]
c = cauchy(V_eta, mpf('0.02'), 40, 6)
l2 = log(2); z2 = pi**2/6; z3 = zeta(3)
pr = {1: -4*l2, 2: 2*l2**2+4*l2, 3: pi**2/6-4*l2-4*l2**2-2*l2**3/3,
      4: z3-pi**2/3-pi**2*l2/6+l2**4/6+2*l2**3+4*l2+6*l2**2,
      5: -8*l2**2-3*z3-4*l2-4*l2**3-l2*z3-2*l2**4/3-l2**5/30+pi**4/720+pi**2*l2**2/12+pi**2*l2/2+pi**2/2,
      6: zeta(5)/2+6*z3+4*l2*z3+l2**2*z3/2-pi**2*z3/12+l2**6/180+l2**5/6+5*l2**4/3+20*l2**3/3+10*l2**2+4*l2-pi**2*l2-2*pi**2/3-29*pi**4/1440-pi**2*l2**2/3-pi**4*l2/720-pi**2*l2**3/36}
print("c_0 =", nstr(c[0], 25))
for k in range(1, 7):
    print("c_%d = %s   notes' exact form diff %s" % (k, nstr(c[k], 25), nstr(c[k]-pr[k], 3)))
# top-weight check: weight-n part of c_n = 4(-ln2)^n/n!  (c_4: ln^4 2/6 = 4 ln^4 2/24 ok)
# circle end: delta = 1/p - 1/2, p = 2/(1+2 delta)
a = cauchy(lambda d: V(2/(1+2*d)), mpf('0.15'), 96, 9)
pa = {2: 3*pi*(pi**2-8)/4, 4: pi*(15*pi**4-352*pi**2+1664)/192,
      6: pi*(9699*pi**6-96280*pi**4+13440*pi**2-39936)/11520,
      8: pi*(-418761675*pi**8+4482752384*pi**6-3451302400*pi**4-63651840*pi**2+138117120)/77414400}
for k in range(1, 10):
    extra = ("  notes' form diff %s" % nstr(a[k]-pa[k], 3)) if k in pa else ""
    print("a_%d = %s%s" % (k, nstr(a[k], 25), extra))
