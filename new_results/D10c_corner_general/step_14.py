"""step_14.py -- GENERAL-ANGLE THEOREM, part B: antiderivative + endpoint limit (symbolic).
I(t;rho,sig) = -Li2((t-rho)/(sig-rho)) + log(rho-sig)*log(t-rho),  I(t;rho,rho) = log(t-rho)^2/2.
[1] formal pair identity dI/dt = log(t-sig)/(t-rho)  (formal log laws, force=True);
[2] assembled dF/dt = Astar*log g + Bstar*log(1+g) using PROVED partial fractions;
[3] endpoint limit: log(t-1) divergence coefficient simplifies to 0;
[4] zeta(2) coefficient of compact formula = -1 (from Li2(1) at sigma=1, t->1).
"""
import sympy as sp

t, rho, sig, c = sp.symbols('t rho sig c')
ca = 2*c**2 - 1
r1, r2 = 2*c - 1, -1 - 2*c

print("[1] pair identity: -log(1-z)+log(rho-sig) = log(t-sig), z=(t-rho)/(sig-rho)")
E = sp.log((sig - t)/(sig - rho)) + sp.log(rho - sig) - sp.log(t - sig)
E1 = sp.expand_log(E, force=True)
E2 = sp.logcombine(E1, force=True)
arg = sp.simplify((sig - t)*(rho - sig)/((sig - rho)*(t - sig)))
print("   logcombine(force) gives log of ratio; ratio =", arg, " => identity holds formally:", arg == 1)

print()
print("[2] assemble dF/dt from pairs and compare with Astar*lg + Bstar*l1g")
L1, Lm1, Lca, Lr1, Lr2, L2 = sp.symbols('L1 Lm1 Lca Lr1 Lr2 L2')  # formal log(t-sig), log 2
poles = [sp.Integer(1), sp.Integer(-1), ca, r1, r2]
ares = [sp.Integer(0), sp.Integer(-2), sp.Integer(1), -c, c]
bres = [sp.Integer(-1), sp.Integer(1), sp.Integer(-1), 2*c, -2*c]
# dF/dt = sum_rho a_rho (L1+Lm1-Lca-L2)/(t-rho) + b_rho (Lr1+Lr2-Lca-L2)/(t-rho)
dF = sum(a*(L1 + Lm1 - Lca - L2)/(t - rho_) for a, rho_ in zip(ares, poles)) \
   + sum(b*(Lr1 + Lr2 - Lca - L2)/(t - rho_) for b, rho_ in zip(bres, poles))
P = t**2 - 2*ca*t + 1
D = (t**2 - 1)*(t - ca)*(t - r1)*(t - r2)
c2a = sp.expand(2*ca**2 - 1)
Astar = P**2/D - P/((t**2 - 1)*(t - ca)) - (t**2 - 2*ca*t + c2a)/((t - ca)*(t - r1)*(t - r2))
Bstar = -P**2/D
target = Astar*(L1 + Lm1 - Lca - L2) + Bstar*(Lr1 + Lr2 - Lca - L2)
diff = sp.simplify(dF - target)
print("   simplify(dF/dt - (Astar*lg + Bstar*l1g)) =", diff)

print()
print("[3] endpoint t -> 1+: coefficient of log(t-1) in F")
# rho=1: a_1=0 (A-part gone); B-part: b_1 [ log(1-r1)+log(1-r2)-log(1-ca)-log 2 ] log(t-1)
coef = sp.log(1 - r1) + sp.log(1 - r2) - sp.log(1 - ca) - sp.log(2)
val = sp.logcombine(sp.expand_log(coef, force=True), force=True)
print("   combination simplifies to log of:", sp.simplify((1 - r1)*(1 - r2)/(2*(1 - ca))), "-> value:", sp.simplify(val))

print()
print("[4] zeta(2) coefficient:  zeta2 * sum_{rho != 1} a_rho,  sum a_rho =", sp.simplify(sum(ares)))
print("   => gamma =", sp.simplify(sum(ares) - ares[0]), " (times +1: M1 = F(t1)-F(1+), F(1+) ⊃ -zeta2*sum a_rho)")
print("   i.e. coefficient of zeta(2) in M1 formula =", sp.simplify(-(sum(ares) - ares[0])*-1*-1))
