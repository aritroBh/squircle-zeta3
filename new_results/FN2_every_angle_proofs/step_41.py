import sympy as sp, time
t0 = time.time()
a, A, B, D4 = sp.symbols('a A B D4')
c, z = sp.symbols('c z')

# theta_ij from previous step (verified vs P3): hard-code the needed ones
th = {
 (1,0): 1/a,
 (2,0): (A - a**2)/(2*a**3),
 (1,1): -A/a**3,
 (3,0): (3*A**2 - 3*A*a**2 - B*a + a**4)/(6*a**5),
 (2,1): (-4*A**2 + 3*A*a**2 + B*a)/(2*a**5),
 (1,2): A**2/a**5,
 (4,0): -(a**6 - 7*a**4*A - 6*a**3*B + 18*a**2*A**2 - a**2*D4 + 10*a*A*B - 15*A**3)/(24*a**7),
 (3,1): -(7*a**4*A + 6*a**3*B - 24*a**2*A**2 + a**2*D4 - 13*a*A*B + 23*A**3)/(6*a**7),
 (2,2): -3*A*(2*a**2*A + a*B - 3*A**2)/(2*a**7),
 (1,3): -A**3/a**7,
}
# check theta_{1j} = a^{-1}(-A/a^2)^j
print("theta13 formula:", sp.simplify(th[(1,3)] - (1/a)*(-A/a**2)**3) == 0)

d2 = 2*z*(1-z); d3 = -4*c*z*(1-z); d4 = 3*(1-c**2)*(5*c**2-1)/8
# T4: l=2: (d2/a)[6 th30 - 2 th21 + th12 + 2 th20 - th11 + th10]
#     l=3: (d3/a^2)[-th20/4 + th11/4 - th10/2]
#     l=4: th10 d4/(9 a^3)
T4 = ( d2/a*(6*th[(3,0)] - 2*th[(2,1)] + th[(1,2)] + 2*th[(2,0)] - th[(1,1)] + th[(1,0)])
     + d3/a**2*(-th[(2,0)]/4 + th[(1,1)]/4 - th[(1,0)]/2)
     + th[(1,0)]*d4/(9*a**3) )
T4 = sp.expand(T4)
T4_P3 = (12*a**4 - 48*a**2*A - a**2 - 24*a*B + 96*A**2)/(24*a**6)
print("T4 at c=0,z=1/2 vs P3:", sp.simplify(T4.subs({c:0, z:sp.Rational(1,2)}) - T4_P3) == 0)
# Also T2, T3 recheck with same code path
T2 = d2/a*th[(1,0)]
T3 = ( d2/a*(-2*th[(2,0)] + th[(1,1)] - th[(1,0)]) + d3/a**2*(th[(1,0)]/4) )
print("T2:", sp.simplify(T2), "\nT3:", sp.simplify(T3))
print("T4 (general c,z):", sp.factor(T4))
print("elapsed", time.time()-t0)
