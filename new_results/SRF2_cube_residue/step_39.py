import mpmath as mp
mp.mp.dps = 100
# Closed-form proof of I = int_0^1 (sqrt(1+a^4)-1)/a^2 da, step by step, all numerically verified.

# (i) I by direct quadrature (stable form a^2/(sqrt(1+a^4)+1))
I_dir = mp.quad(lambda a: a*a/(mp.sqrt(1+a**4)+1), [0, 1])
print("I direct          =", mp.nstr(I_dir, 40))

# (i') integration by parts: I = 1 - sqrt(2) + 2K, K = int_0^1 a^2/sqrt(1+a^4) da
K = mp.quad(lambda a: a*a/mp.sqrt(1+a**4), [0, 1])
print("I via IBP         =", mp.nstr(1-mp.sqrt(2)+2*K, 40))

# (ii) L' = int_0^1 t^2/sqrt(1-t^4) dt = (1/4) B(3/4,1/2) = sqrt(2) pi^{3/2}/Gamma(1/4)^2
Lp_dir = mp.quad(lambda t: t*t/mp.sqrt(1-t**4), [0, 1])
Lp_beta = mp.mpf('0.25')*mp.beta(mp.mpf('0.75'), mp.mpf('0.5'))
Lp_closed = mp.sqrt(2)*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2
print("L' direct         =", mp.nstr(Lp_dir, 40))
print("L' = (1/4)B(3/4,1/2) =", mp.nstr(Lp_beta, 40))
print("L' closed         =", mp.nstr(Lp_closed, 40))

# (iii) substitution t = sqrt(2) a/sqrt(1+a^4):  L' = 2 sqrt(2) int_0^1 a^2/(1+a^4)^{3/2} da
Lp_sub = 2*mp.sqrt(2)*mp.quad(lambda a: a*a/(1+a**4)**mp.mpf('1.5'), [0, 1])
print("L' via subst      =", mp.nstr(Lp_sub, 40))

# (iv) int_0^1 a^2/(1+a^4)^{3/2} da = 1/(2 sqrt(2)) - K/2   [derivative identity]
lhs = mp.quad(lambda a: a*a/(1+a**4)**mp.mpf('1.5'), [0, 1])
rhs = 1/(2*mp.sqrt(2)) - K/2
print("IBP2 lhs          =", mp.nstr(lhs, 40))
print("IBP2 rhs          =", mp.nstr(rhs, 40))

# (v) K = (1 - L')/sqrt(2)
print("K vs (1-L')/sqrt2 =", mp.nstr(K, 40), " ", mp.nstr((1-Lp_closed)/mp.sqrt(2), 40))

# (vi) final: I = 1 - 2 pi^{3/2}/Gamma(1/4)^2
I_closed = 1 - 2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2
print("I closed          =", mp.nstr(I_closed, 40))
print("R = I - 2         =", mp.nstr(I_closed-2, 40))
print("planner R         =", mp.nstr(-1-2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2, 40))
