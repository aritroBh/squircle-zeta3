import mpmath as mp
mp.mp.dps = 50
# Double series Li_{m,n}(x,y) = sum_{p>q>=1} x^p y^q / (p^m q^n)
def li_mn(m, n, x, y):
    return mp.fsum(x**p * y**q / (p**m * q**n) for p in range(1, 400) for q in range(1, p))
# Iterated integral G_0(0^{m-1}, 1/x, 0^{n-1}, 1/(x y); 1), base point 0 -> divergent at 0; use standard:
# G_0(word;1) with word starting with 0's: this is -log divergence... Actually MPL iterated integrals
# I(0; w; 1) converge iff last letter != 0. Word ends in 1/(xy) != 0 so convergent.
# G_c(w; e) = integral_c^e du/(u - w_1) G_c(w[1:]; u), G_c((); u) = 1. With c=0, first letters 0 diverge?
# No: the standard definition handles it: iterated integral int_0^1 dt1/(t1-a1) ... with a1=0: log 0 divergence!
# BUT the MPL definition uses the convention with REVERSED order: Li_{m,n}(x,y) = (-1)^k I(0; a; 1) where
# the INNERMOST integration is at 0. Let me just test both sign/word conventions numerically via limits.
def G(word, e, c=mp.mpf(0)):
    # integral from c to e, first letter outermost; regularize at c by eps-limit
    if not word:
        return mp.mpf(1)
    a = word[0]
    # int_c^e du/(u-a) G(word[1:], u)
    # near u=c: if G(word[1:],c) finite and a != c, fine; if a == c == 0, divergence unless G(word[1:],c)=0
    f = lambda u: G(word[1:], u, c)/(u-a) if abs(u-a) > mp.mpf('1e-40') else mp.mpf(0)
    return mp.quad(f, [c, e], maxdegree=8)

x, y = mp.mpf(1)/3, mp.mpf(1)/2
for (m, n) in [(3, 2), (2, 1)]:
    word = (mp.mpf(0),)*(m-1) + (1/x,) + (mp.mpf(0),)*(n-1) + (1/(x*y),)
    # check innermost: G(();u)=1; G((1/(xy),);u) = int_0^u dt/(t - 1/(xy)) = log(1 - xy*u)  [finite]
    # next letters 0: int dt/t log(...) finite. Outermost 0's fine. So no divergence.
    val = G(word, mp.mpf(1))
    ser = li_mn(m, n, x, y)
    print(f"Li_{{{m},{n}}}({x},{y}): series={mp.nstr(ser,40)}")
    print(f"   G_0(word;1) = {mp.nstr(val,40)}, ratio = {mp.nstr(val/ser, 20)}")