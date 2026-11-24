
import mpmath as mp
mp.mp.dps = 50

def prepend_a(s, a, N):
    """J((a,rest); t) from series s of J(rest; t): t_new[n+1] = -acc_n/(a^{n+1}(n+1))."""
    t = [mp.mpf(0)] * (N + 1)
    aP = mp.mpf(1)
    acc = mp.mpf(0)
    for n in range(1, N):
        aP *= a
        acc += s[n] * aP
        t[n + 1] = -acc / (aP * a * (n + 1))
    return t

# test: J((2,0,3); 1) = ∫₀¹ du/(u-2) * (-Li2(u/3))
N = 800
a, b = mp.mpf(2), mp.mpf(3)
s = [mp.mpf(0)] * (N + 1)
for m in range(1, N + 1):
    s[m] = -(mp.mpf(1)/b)**m / m
s = [mp.mpf(0)] + [s[m]/m for m in range(1, N+1)]   # prepend 0
t = prepend_a(s, a, N)
J_series = mp.fsum(t[1:])
J_direct = mp.quad(lambda u: -mp.polylog(2, u/3)/(u - 2), [0, 1])
print("J((2,0,3);1) series:", mp.nstr(J_series, 20))
print("J((2,0,3);1) direct:", mp.nstr(J_direct, 20))
print("diff:", mp.nstr(abs(J_series - J_direct), 3))
