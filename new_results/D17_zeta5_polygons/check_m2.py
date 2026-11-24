
import sympy as sp, sys, time
import mpmath as mp
mp.mp.dps = 80
sys.path.insert(0, "research/08_theory/D17_zeta5_polygons")
import sym_D5

# Build reduced D5, extract cross G's, lambdify
D, syms = sym_D5.build_D(5, reduced=True)
L0,L1 = syms['L0'], syms['L1']
s0,s1,c0,c1,S0,Sc = syms['s0'],syms['s1'],syms['c0'],syms['c1'],syms['S0'],syms['Sc']
D5 = sp.expand(D[5])
Gs = {}
for (a,b) in [(4,1),(3,2),(2,3),(1,4)]:
    Gs[(a,b)] = sp.expand(D5.coeff(L0,a).coeff(L1,b))
Gsum = sp.expand(Gs[(4,1)]+Gs[(3,2)]+Gs[(2,3)]+Gs[(1,4)])
f = sp.lambdify((s0,s1,c0,c1,S0,Sc), Gsum, 'mpmath')

# m=2: on [0, pi/2], s0=sin, s1=cos, c0=cos, c1=-sin, S0=sin+cos, Sc=cos-sin
def vals_m2(th):
    s0v=mp.sin(th); s1v=mp.cos(th); c0v=mp.cos(th); c1v=-mp.sin(th)
    return (s0v,s1v,c0v,c1v, s0v+s1v, c0v+c1v)
A = mp.pi/2
J5 = mp.quad(lambda th: f(*vals_m2(th)), [0, A], method='tanh-sinh')
print("dps:", mp.mp.dps)
print("J5(m=2) = int_0^{pi/2} (G41+G32+G23+G14) =", mp.nstr(J5, 40))
target = -mp.sqrt(2)/128
print("target -sqrt(2)/128 =", mp.nstr(target, 40))
print("RESIDUAL =", mp.nstr(abs(J5-target), 5))
# Then C_2^(5) = 5760*2*J5 = 11520 J5 ; target -90 sqrt(2)
C2 = 11520*J5
print("C_2^(5) = 11520*J5 =", mp.nstr(C2, 30), " target -90sqrt2 =", mp.nstr(-90*mp.sqrt(2),30))
print("C2 residual:", mp.nstr(abs(C2+90*mp.sqrt(2)),5))
