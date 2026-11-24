# ---- relabel zeta symbols, build a_n in d, reduce even zetas to pi ----
old = list(z)                                    # old symbols z2a16,... (order = zeta(2),zeta(3),...)
Z = sp.symbols(' '.join(f'Z{n}' for n in range(2, 2+len(old))))
remap = dict(zip(old, Z))
def zet(n): return Z[n-2]

Fc_c = [sp.expand(f.subs(remap)) for f in Fc]
# P/sqrt2 = 4(1+eta) F(eta); eta = -d/(1+d)
Ac_eta = [sp.expand(4*(Fc_c[m] + (Fc_c[m-1] if m >= 1 else 0))) for m in range(0, 11)]
# substitute eta(d) = sum_{j>=1} (-1)^j d^j, keep to d^10
eta_pow = [sp.Integer(1)]
for m in range(1, 11):
    eta_pow.append(sp.expand(-T(eta_pow[-1], d, 10)))   # eta^m = -d * eta^{m-1} truncated (since eta = -d*(1 - d + ...))
# careful: eta = -d + d^2 - ... ; multiply by geometric: eta^m = (-d)^m * (1-d+...)^m. Build directly:
eta_ser = sum((-d)**j for j in range(1, 12))
eta_pow = [sp.Integer(1)]
for m in range(1, 11):
    eta_pow.append(T(eta_pow[-1]*eta_ser, d, 10))
Ac = [sp.expand(sum(Ac_eta[m]*eta_pow[m] for m in range(0, 11))).coeff(d, n) for n in range(0, 11)]

pi_red = {zet(2): sp.pi**2/6, zet(4): sp.pi**4/90, zet(6): sp.pi**6/945, zet(8): sp.pi**8/9450, zet(10): sp.pi**10/93555}
print("a_0 =", sp.expand(Ac[0].subs(pi_red)))
print("a_1 =", sp.expand(Ac[1].subs(pi_red)))
print("a_2 =", sp.expand(Ac[2].subs(pi_red)), "   [planner: pi^2/6]")
print("a_3 =", sp.expand(Ac[3].subs(pi_red)), "   [planner: zeta(3) - pi^2/3]")
for n in range(4, 9):
    print(f"a_{n} =", sp.expand(Ac[n].subs(pi_red)))
print()
for n in range(9, 11):
    print(f"a_{n} =", sp.expand(Ac[n].subs(pi_red)))
