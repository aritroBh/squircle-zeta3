# SM5: numerical check (dps 60, quad) of every integral used in the exact derivation.
import mpmath as mp
mp.mp.dps = 60
z2, z3, z5 = mp.zeta(2), mp.zeta(3), mp.zeta(5)
g  = lambda s: mp.log(1 + mp.exp(-2*s))
gp = lambda s: -2/(1 + mp.exp(2*s))
def gam(a, b):
    return mp.quad(lambda s: s**a*g(s)**b, [0, 1, 4, 16, 64, mp.inf])/(mp.factorial(a)*mp.factorial(b))
G = {(a, 4-a): gam(a, 4-a) for a in range(0, 4)}
print("gamma_ab (a+b=4):"); [print(" ", k, mp.nstr(v, 45)) for k, v in G.items()]
chk = {
 "gam31 = 15 z5/256": G[(3,1)] - mp.mpf(15)/256*z5,
 "gam22 = -3 gam31 + z2z3/16 + z5/16": G[(2,2)] - (-3*G[(3,1)] + z2*z3/16 + z5/16),
 "gam04 = -gam13 + 2 gam31 - z2z3/8 + z5/8": G[(0,4)] - (-G[(1,3)] + 2*G[(3,1)] - z2*z3/8 + z5/8),
}
for k, v in chk.items(): print("residual", k, mp.nstr(v, 3))
# soft-max: 8 * int_0^inf Phi4, Phi4 from sm5_derive_output.txt, g3 = s^3/(3(1+e^{2s}))
g3 = lambda s: s**3/(3*(1 + mp.exp(2*s)))
Phi4 = lambda s: g(s)**4/24 + g(s)**3*s/6 - g(s)**2*s**2/4 + g(s)*g3(s) - g(s)*s**3/6 + g3(s)*s
I4 = mp.quad(Phi4, [0, 1, 4, 16, 64, mp.inf])
print("softmax: 8*int Phi4 =", mp.nstr(8*I4, 50))
print("closed form 107/32 z5 - z2 z3 =", mp.nstr(mp.mpf(107)/32*z5 - z2*z3, 50))
print("residual:", mp.nstr(8*I4 - (mp.mpf(107)/32*z5 - z2*z3), 3))
print("Richardson c5 (sm5_c5_numeric.txt) minus closed form:", mp.nstr(mp.mpf(open('sm5_c5_numeric.txt').read().strip()) - (mp.mpf(107)/32*z5 - z2*z3), 3))
# squircle: Phi4 with H=1/sqrt2, g3 = -4 s^3/(3(1+e^{2s}))
r2 = mp.sqrt(2)
g3q = lambda s: -4*s**3/(3*(1 + mp.exp(2*s)))
def Phi4q(s):
    G1 = g(s); G2 = G1**2/2; G3 = G1**3/6 + g3q(s); G4 = G1**4/24 + G1*g3q(s)
    return (G4 + s*G3 - s**2/2*G2 - s**3/6*G1)/r2
Iq = mp.quad(Phi4q, [0, 1, 4, 16, 64, mp.inf])
print("squircle: 8*int Phi4 =", mp.nstr(8*Iq, 50))
print("-sqrt2(3/4 z5 + 7/4 z2 z3) =", mp.nstr(-r2*(mp.mpf(3)/4*z5 + mp.mpf(7)/4*z2*z3), 50))
print("residual:", mp.nstr(8*Iq + r2*(mp.mpf(3)/4*z5 + mp.mpf(7)/4*z2*z3), 3))
