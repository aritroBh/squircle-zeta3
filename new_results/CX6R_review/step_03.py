import mpmath as mp
mp.mp.dps = 50
# CX6R step_03: test proof.md eq (3) (logit form of M(eps)) against the
# original definition M = -F/2 (paper eq:F), for several real eps.
# Original: M(eps) = -int_0^{1/2} (a+b-sqrt(a^2+b^2)) dt,
#   a = t^{eps-1}, b = (1-t)^{eps-1}, sqrt = a*sqrt(1+(b/a)^2) (branch of lem:holo).
def M_orig(eps):
    def f(t):
        a = t**(eps-1); b = (1-t)**(eps-1)
        return a + b - a*mp.sqrt(1+(b/a)**2)
    return -mp.quad(f, [mp.mpf('1e-40'), mp.mpf('0.25'), mp.mpf('0.5')])
def B(z):
    return 1 + mp.exp(-z) - mp.sqrt(1+mp.exp(-2*z))
def C(z):
    return mp.exp(z)*B(z)
def M_form1(eps):  # -int e^{-eps x}(1+e^{-x})^{-1-eps} B((1-eps)x) dx
    return -mp.quad(lambda x: mp.exp(-eps*x)*(1+mp.exp(-x))**(-1-eps)*B((1-eps)*x),
                    [0, 1, 10, 100, mp.inf])
def M_form2(eps):  # -int e^{-x}(1+e^{-x})^{-1-eps} C((1-eps)x) dx
    return -mp.quad(lambda x: mp.exp(-x)*(1+mp.exp(-x))**(-1-eps)*C((1-eps)*x),
                    [0, 1, 10, 100, mp.inf])
print("CX6R step_03: eq (3) vs original definition, 50 dps")
worst = mp.mpf(0)
for eps in [mp.mpf('0'), mp.mpf('0.1'), mp.mpf('-0.2'), mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.9')]:
    m0 = M_orig(eps); m1 = M_form1(eps); m2 = M_form2(eps)
    d1 = abs(m0-m1); d2 = abs(m0-m2)
    worst = max(worst, d1, d2)
    print(f"eps={mp.nstr(eps,6)}  M={mp.nstr(m0,30)}")
    print(f"   |orig-form1|={mp.nstr(d1,6)}  |orig-form2|={mp.nstr(d2,6)}")
print("worst abs discrepancy:", mp.nstr(worst, 6))
