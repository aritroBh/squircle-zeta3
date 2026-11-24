
import mpmath as mp, json
mp.mp.dps = 50

# V(eps) analytic continuation to Re eps > -1 via subtractions:
# V = V_sub + Pa + Pc
# V_sub = int_0^1 int_0^1 a^{2e-1} c^{e-1} [g - g(a,0) - a A1(c) - c B1(a)] da dc
# A1(c) = -(1+2e)[(1+c)^{-2-2e} - 1]
# B1(a) = -(1+2e) a (1+a)^{-2-2e} sqrt(1+a^{2-2e})
# Pa = -(1+2e)/(2e+1) * int_0^1 c^{e-1}[(1+c)^{-2-2e}-1] dc
# Pc = -(1+2e)/(e+1) * int_0^1 a^{2e}(1+a)^{-2-2e} sqrt(1+a^{2-2e}) da

def Vcont(e):
    e = mp.mpf(e)
    def g(a, c):
        N = 1 + a + a*c
        return N**(-1-2*e)*mp.sqrt(1 + a**(2-2*e) + (a*c)**(2-2*e))
    def integrand(a, c):
        A1 = -(1+2*e)*((1+c)**(-2-2*e) - 1)
        B1 = -(1+2*e)*a*(1+a)**(-2-2*e)*mp.sqrt(1+a**(2-2*e))
        return a**(2*e-1)*c**(e-1)*(g(a,c) - g(a,mp.mpf(0)) - a*A1 - c*B1)
    Vsub = mp.quad(integrand, [0,1], [0,1], method='tanh-sinh')
    Ja = mp.quad(lambda c: c**(e-1)*((1+c)**(-2-2*e)-1), [0,1])
    Pa = -(1+2*e)/(2*e+1)*Ja
    Ic = mp.quad(lambda a: a**(2*e)*(1+a)**(-2-2*e)*mp.sqrt(1+a**(2-2*e)), [0,1])
    Pc = -(1+2*e)/(e+1)*Ic
    return Vsub + Pa + Pc

# check pole: (e+1)*V(e) as e -> -1
print(" e        V(e)            (e+1)V(e)")
for e in [-0.7, -0.8, -0.85, -0.9, -0.95]:
    v = Vcont(e)
    print(f"{e: .2f}  {mp.nstr(v,10)}  {mp.nstr((e+1)*v,10)}")
