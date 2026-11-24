# step_50: control: eq:chi3 (proved in the paper) vs the single-valued L3 framework
import mpmath as mp
mp.mp.dps = 120
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
z3 = mp.zeta(3)
# L(chi8,3): chi8(1)=chi8(7)=1, chi8(3)=chi8(5)=-1
def chi8(n):
    r = n % 8
    return 0 if r%2==0 else (1 if r in (1,7) else -1)
L8 = mp.nsum(lambda n: chi8(int(n))/mp.mpf(n)**3, [1, mp.inf])
print("L(chi8,3) =", L8, " paper: 0.95838045456...")
b = mp.sqrt(2)-1
L = mp.log(1+mp.sqrt(2))
# eq:chi3: 48 chi_3(b) = 32 sqrt2 L(chi8,3) - 3 pi^2 L + 4 L^3
chi3b = (Li3(b)-Li3(-b))/2
lhs = 48*chi3b
rhs = 32*mp.sqrt(2)*L8 - 3*mp.pi**2*L + 4*L**3
print("eq:chi3 check: lhs-rhs =", lhs-rhs)
# now the L3-class value: D8 = L3(b)-L3(-b); is it q*L(chi8,3)?
D8 = L3(b)-L3(-b)
print("D8 =", D8, " D8/L8 =", D8/L8)
import sympy as sp
print("D8/L8 nsimplify:", sp.nsimplify(str(D8/L8), [sp.sqrt(2)]))
# predicted: {b}-{-b} = 2 beta8 {chi8}; chi3(b) = (2sqrt2/3)L8 - pi^2 L/16 + L^3/12
# D8 = 2 chi3(b) - logb (Li2(b)-Li2(-b)) - (1/3) log^2 b log((1-b)/(1+b))
cor = -mp.log(b)*(Li2(b)-Li2(-b)) - mp.mpf(1)/3*mp.log(b)**2*mp.log((1-b)/(1+b))
print("correction terms:", cor)
print("D8 - (4 sqrt2/3) L8 =", D8 - 4*mp.sqrt(2)/3*L8)
