
import mpmath as mp
mp.mp.dps = 105
r5 = mp.sqrt(5); phi = (1+r5)/2
s = mp.sin(mp.pi/10); c = mp.cos(mp.pi/10)
C5_5 = 2880*5*s**4*c*(65*s**6 - 63*s**4 + 17*s**2 - 1)
C5_23 = -960*5*c*s**4*(45*s**6 - 27*s**4 + 6*s**2 - 1)
L = mp.sqrt((25-11*r5)/10)
print("C5^(5)  =", mp.nstr(C5_5, 40))
print("C5^(2,3)=", mp.nstr(C5_23, 40))
# predicted local-block zeta coefficients in P5/L: divide by L
print("C5^(5)/L  =", mp.nstr(C5_5/L, 40))
print("C5^(2,3)/L=", mp.nstr(C5_23/L, 40))
# are C5^(5)/L and C5^(2,3)/L in Q(sqrt5)?
import sympy as sp
x = sp.nsimplify(C5_5/L, [sp.sqrt(5)])
print("nsimplify C5^(5)/L:", x)
y = sp.nsimplify(C5_23/L, [sp.sqrt(5)])
print("nsimplify C5^(2,3)/L:", y)
