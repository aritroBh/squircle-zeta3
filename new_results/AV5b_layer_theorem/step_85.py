import mpmath as mp
mp.mp.dps = 60
g3 = lambda s: (abs(s)**3/6)*(1-mp.tanh(s))
Ig3 = 2*mp.quad(g3, [0, 30])
print("int_R g3 =", mp.nstr(Ig3, 40))
print("7 zeta(4)/32 =", mp.nstr(7*mp.zeta(4)/32, 40))
print("resid:", mp.nstr(Ig3 - 7*mp.zeta(4)/32, 8))
# c4 for F1 = 4[Ig3 + 2J - int s^2 g1/2], J = zeta(4)/32
J = mp.zeta(4)/32
Is2g1 = mp.quad(lambda s: s**2*mp.log(1+mp.exp(-2*s)), [0, mp.inf])
c4 = 4*(Ig3 + 2*J - Is2g1)
print("c4_F1 predicted =", mp.nstr(c4, 40))
print("measured A4 = 0.270590261081837161503305710993")
print("zeta(4)/4 =", mp.nstr(mp.zeta(4)/4, 40))
