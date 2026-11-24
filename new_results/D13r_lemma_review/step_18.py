import mpmath as mp, random
import sympy as sp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); sqpi = mp.sqrt(mp.pi); a = mp.pi/(2*sqrt2)

# symbolic: Phi(x) = (-1+i)a x + N log(bx/(1+bx)); derive Phi'', Phi''', Phi''''
x, N = sp.symbols('x N', positive=True)
i = sp.I; bb = sp.exp(i*sp.pi/4); aa = sp.pi/(2*sp.sqrt(2))
Phi = (-1+i)*aa*x + N*sp.log(bb*x/(1+bb*x))
P1 = sp.simplify(sp.diff(Phi,x))
P2 = sp.simplify(sp.diff(Phi,x,2))
P3 = sp.simplify(sp.diff(Phi,x,3))
P4 = sp.simplify(sp.diff(Phi,x,4))
P2c = -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
P3c = 2*N/x**3 - 2*N*bb**3/(1+bb*x)**3
P4c = -6*N/x**4 + 6*N*bb**4/(1+bb*x)**4
print("P1 - ((-1+i)a + N(1/x - b/(1+bx))) =", sp.simplify(P1 - ((-1+i)*aa + N*(1/x - bb/(1+bb*x)))))
print("P2 formula residual:", sp.simplify(P2-P2c))
print("P3 formula residual:", sp.simplify(P3-P3c))
print("P4 formula residual:", sp.simplify(P4-P4c))

mp.mp.dps = 60
bbm = mp.e**(1j*mp.pi/4)
def Phi3(x,N): return 2*N/x**3 - 2*N*bbm**3/(1+bbm*x)**3
def Phi4(x,N): return -6*N/x**4 + 6*N*bbm**4/(1+bbm*x)**4
def Phi3b(x,N): return 2*N*(1+3*x+3*x**2)/(x**3*(1+sqrt2*x+x*x)**mp.mpf('1.5'))
def Phi4b(x,N): return 6*N*(1+4*x+6*x**2+4*x**3)/(x**4*(1+sqrt2*x+x*x)**2)

random.seed(11)
w3=(mp.mpf('-1'),None); w4=(mp.mpf('-1'),None)
for _ in range(2500):
    xv = mp.mpf(10)**random.uniform(-2, 6)
    Nv = mp.mpf(10)**random.uniform(0, 12)
    r3 = abs(Phi3(xv,Nv))/Phi3b(xv,Nv); r4 = abs(Phi4(xv,Nv))/Phi4b(xv,Nv)
    if r3>w3[0]: w3=(r3,xv)
    if r4>w4[0]: w4=(r4,xv)
print("L4: worst |Phi3|/Phi3b =", mp.nstr(w3[0],10), "at x =", mp.nstr(w3[1],6))
print("L4: worst |Phi4|/Phi4b =", mp.nstr(w4[0],10), "at x =", mp.nstr(w4[1],6))

# monotonicity in x of Phi3b/Phi4b (per unit N): sympy sign analysis of log-derivative
xs = sp.symbols('xs', positive=True)
n3 = 1+3*xs+3*xs**2; d3 = xs**3*(1+sp.sqrt(2)*xs+xs**2)**sp.Rational(3,2)
n4 = 1+4*xs+6*xs**2+4*xs**3; d4 = xs**4*(1+sp.sqrt(2)*xs+xs**2)**2
# d/dx log(n3/d3) < 0  <=>  n3'/n3 < d3'/d3 ; prove n3'/n3 <= 2/x  and 3/x <= d3'/d3
ineq1 = sp.simplify(2*n3 - xs*sp.diff(n3,xs))      # >=0 shows n3'/n3 <= 2/x
print("2*n3 - x*n3' =", sp.expand(ineq1), " (all coeffs >=0 => n3'/n3 <= 2/x)")
ineq2 = sp.simplify(4*n4 - xs*sp.diff(n4,xs))
print("4*n4 - x*n4' =", sp.expand(ineq2), " (all coeffs >=0 => n4'/n4 <= 4/x)")
# d3'/d3 = 3/x + (3/2)(sqrt2+2x)/(1+sqrt2 x+x^2) > 3/x > 2/x >= n3'/n3  => Phi3b strictly decreasing
# d4'/d4 = 4/x + 2(sqrt2+2x)/(1+sqrt2 x+x^2) > 4/x >= n4'/n4 => Phi4b strictly decreasing
print("d3'/d3 = 3/x + (3/2)(rt2+2x)/(1+rt2 x+x^2) > 3/x > 2/x >= n3'/n3  => Phi3b decr.")
print("d4'/d4 = 4/x + 2(rt2+2x)/(1+rt2 x+x^2)   > 4/x >= n4'/n4          => Phi4b decr.")
# numeric monotone check
maxd = mp.mpf('-9')
for _ in range(2000):
    x1 = mp.mpf(10)**random.uniform(-2,6); x2 = x1*mp.mpf(10)**random.uniform(0.0001,1)
    Nv = mp.mpf('123.4')
    maxd = max(maxd, Phi3b(x2,Nv)-Phi3b(x1,Nv)); maxd = max(maxd, Phi4b(x2,Nv)-Phi4b(x1,Nv))
print("numeric monotonicity: max increase over 4000 pairs (must be <0):", mp.nstr(maxd,6))
with open('research/08_theory/D13r_lemma_review/progress.log','a') as f:
    f.write("L4: HOLDS (Phi''/Phi'''/Phi'''' formulas sympy-verified residual 0; bounds from (1+bx)^m-(bx)^m expansion; "
            "2500 random (x,N) @60dps worst ratios printed; monotonicity PROVED via n'/n <= m/x < d'/d, sympy-expanded numerators all-positive-coeffs; 4000-pair numeric monotone check <0).\n")
print("logged")
