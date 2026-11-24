"""
F1 step 1: the perimeter of the p-squircle as the mean of its support function (Cauchy's formula).

Claim (classical ingredients): for p >= 1 the support function of B_p = {|x|^p+|y|^p <= 1} is the dual norm
    h_p(theta) = (|cos theta|^q + |sin theta|^q)^(1/q),   1/p + 1/q = 1,
and Cauchy's formula gives P(p) = int_0^{2pi} h_p(theta) dtheta = pi * a_0(h_p).
With u = tan(theta) on [0, pi/4] (8 symmetric copies):
    P(p) = 8 int_0^1 (1+u^q)^(1/q) (1+u^2)^(-3/2) du,
and with u = exp(-tau*s), tau = 1/q = 1 - 1/p:
    P(p) = 4 sqrt2 + 8 tau int_0^inf [(1+e^-s)^tau - 1] e^(-tau s) (1+e^(-2 tau s))^(-3/2) ds.

We compare against an independent arclength formula (the note's Section 2):
    P/4 = 2 - eps*F(eps),  F = int_0^1 (a+b-sqrt(a^2+b^2)) dt,  a=t^(eps-1), b=(1-t)^(eps-1), eps=1/p.
"""
import mpmath as mp

mp.mp.dps = 50

def P_arclength(p):
    eps = mp.mpf(1) / p
    def phi(t):
        a = t ** (eps - 1); b = (1 - t) ** (eps - 1)
        # stable form: a+b-sqrt(a^2+b^2) = 2ab/(a+b+sqrt(a^2+b^2))
        return 2 * a * b / (a + b + mp.sqrt(a * a + b * b))
    F = 2 * mp.quad(phi, [0, mp.mpf(1) / 4, mp.mpf(1) / 2])
    return 4 * (2 - eps * F)

def P_support_u(p):
    q = p / (p - 1)
    f = lambda u: (1 + u ** q) ** (1 / q) * (1 + u * u) ** (-mp.mpf(3) / 2)
    pts = [0, mp.mpf(1) / 2, 1 - 4 / q, 1 - 1 / q, 1] if q > 8 else [0, mp.mpf(1) / 2, 1]
    return 8 * mp.quad(f, pts)

def P_support_s(p):
    tau = 1 - mp.mpf(1) / p
    g = lambda s: ((1 + mp.e ** (-s)) ** tau - 1) * mp.e ** (-tau * s) * (1 + mp.e ** (-2 * tau * s)) ** (-mp.mpf(3) / 2)
    return 4 * mp.sqrt(2) + 8 * tau * mp.quad(g, [0, 1, 5, 20, 80, mp.inf])

def P_support_theta(p):
    q = p / (p - 1)
    h = lambda th: (abs(mp.cos(th)) ** q + abs(mp.sin(th)) ** q) ** (1 / q)
    return 8 * mp.quad(h, [0, mp.pi / 8, mp.pi / 4])

print("p, P_arclength, P_support(u), P_support(s), P_support(theta), max |diff|")
for p in [mp.mpf('1.05'), mp.mpf('1.25'), mp.mpf('1.5'), mp.mpf(2), mp.mpf(3), mp.mpf(7), mp.mpf(30)]:
    A = P_arclength(p); B = P_support_u(p); C = P_support_s(p); D = P_support_theta(p)
    d = max(abs(A - B), abs(A - C), abs(A - D))
    print(mp.nstr(p, 5), mp.nstr(A, 30), mp.nstr(B, 12), mp.nstr(C, 12), mp.nstr(D, 12), mp.nstr(d, 3))

print("P(2) - 2pi =", mp.nstr(P_support_u(mp.mpf(2)) - 2 * mp.pi, 3))
# endpoint values
print("q->inf (p=1): 4sqrt2 =", mp.nstr(4 * mp.sqrt(2), 20), " ; q=1 (p=inf): 8")

# first-order term at p = infinity: dP/dq at q = 1 should be 8*M0 = 8(log2 - sqrt2*log(1+sqrt2))
M0 = mp.log(2) - mp.sqrt(2) * mp.log(1 + mp.sqrt(2))
def P_of_q(q):
    h = lambda th: (abs(mp.cos(th)) ** q + abs(mp.sin(th)) ** q) ** (1 / q)
    return 8 * mp.quad(h, [0, mp.pi / 8, mp.pi / 4])
dPdq = mp.diff(P_of_q, 1, direction=1)
print("dP/dq at q=1 (one-sided):", mp.nstr(dPdq, 30), "  8*M0 =", mp.nstr(8 * M0, 30), " diff", mp.nstr(dPdq - 8 * M0, 3))
# exact three-line evaluation of that derivative:
# P'(q=1) = 4 int_0^{pi/2} [c log c + s log s - (c+s)log(c+s)] = 4[2(log2-1) - (2sqrt2 L - 2)] = 8(log2 - sqrt2 L)
L = mp.log(1 + mp.sqrt(2))
I1 = mp.quad(lambda t: mp.cos(t) * mp.log(mp.cos(t)), [0, mp.pi / 2])
I2 = mp.quad(lambda t: (mp.cos(t) + mp.sin(t)) * mp.log(mp.cos(t) + mp.sin(t)), [0, mp.pi / 2])
print("int c log c =", mp.nstr(I1, 25), " vs log2-1 =", mp.nstr(mp.log(2) - 1, 25))
print("int (c+s)log(c+s) =", mp.nstr(I2, 25), " vs 2sqrt2 L - 2 =", mp.nstr(2 * mp.sqrt(2) * L - 2, 25))

# derivative at p = 2:  dP/dp|_{p=2} = pi(log 2 - 1/2)
dPdp2 = mp.diff(P_support_u, mp.mpf(2))
print("dP/dp at p=2:", mp.nstr(dPdp2, 30), " pi(log2-1/2) =", mp.nstr(mp.pi * (mp.log(2) - mp.mpf(1) / 2), 30),
      " diff", mp.nstr(dPdp2 - mp.pi * (mp.log(2) - mp.mpf(1) / 2), 3))
