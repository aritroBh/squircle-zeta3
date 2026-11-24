"""Referee's independent computation of I_n = [eta^n] I(eta), I = P/(4 r (1+eta)).

Method (independent of the author's Lemma 2 / u-representation):
  P = 4 (1+eta) int_0^1 sqrt(A^2 t^{2 eta} + B^2 (1-t)^{2 eta}) dt   (arc x=A t^{1/p}, y=B(1-t)^{1/p})
so I(eta) = (1/r) int_0^1 sqrt(A^2 e^{2 eta a} + B^2 e^{2 eta b}) dt,  a=log t, b=log(1-t).
At each quadrature node the integrand is expanded as a power series in eta (series square root),
and the coefficients are integrated by a hand-written tanh-sinh rule, t = 1/(1+e^{-v}), v = pi sinh x.
Convergence is checked by halving the step.
"""
import sys
from mpmath import mp, mpf, sqrt, log, log1p, exp, sinh, cosh, pi, factorial

def sqrt_series(c, N):
    s = [sqrt(c[0])]
    for n in range(1, N + 1):
        acc = c[n] - sum(s[k] * s[n - k] for k in range(1, n))
        s.append(acc / (2 * s[0]))
    return s

def In_direct(A, B, N, h, X=mpf(5.6)):
    A, B = mpf(A), mpf(B)
    r = sqrt(A * A + B * B)
    tot = [mpf(0)] * (N + 1)
    k = -int(X / h)
    while k * h <= X:
        x = k * h
        v = pi * sinh(x)
        a = -log1p(exp(-v))          # log t
        b = -log1p(exp(v))           # log(1-t)
        w = exp(a + b) * pi * cosh(x) * h   # t(1-t) dv/dx h
        c = []
        pa, pb, f = mpf(1), mpf(1), mpf(1)
        for j in range(N + 1):
            if j:
                pa *= 2 * a; pb *= 2 * b; f *= j
            c.append((A * A * pa + B * B * pb) / f)
        s = sqrt_series(c, N)
        for j in range(N + 1):
            tot[j] += s[j] * w
        k += 1
    return [t / r for t in tot]

if __name__ == "__main__":
    mp.dps = 50
    N = 40
    import pickle, os
    out = {}
    for (A, B) in [(1, 1), (1, 3), (2, 5)]:
        I1 = In_direct(A, B, N, mpf(1) / 32)
        I2 = In_direct(A, B, N, mpf(1) / 64)
        err = max(abs(I1[n] - I2[n]) / abs(I2[n]) for n in range(N + 1))
        print(f"(A,B)=({A},{B}) step-halving max rel diff over n<=40: {mp.nstr(err, 3)}", flush=True)
        out[(A, B)] = [mp.nstr(x, 45) for x in I2]
    with open(os.path.join(os.path.dirname(__file__), "In_direct.pkl"), "wb") as fh:
        pickle.dump(out, fh)
