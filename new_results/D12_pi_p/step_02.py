# step_02.py -- TASK 2 setup: exact split pi_p = 4*2^(-1/p) + R(p), and J-integral evaluations
import mpmath as mp
mp.mp.dps = 80
pi = mp.pi; z3 = mp.zeta(3)

def pi_cosh(p):
    p = mp.mpf(p); h = 1/p
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f, [0, mp.inf])

def split(p):  # 4*2^(-h) + R(h), R via s-integral (Lemma 2 of writeup)
    p = mp.mpf(p); h = 1/p
    A = 4*mp.power(2,-h)
    w = lambda s: mp.exp(-s*h*h/(1-h))*(1+mp.exp(-s*h/(1-h)))**(-1-h)
    R = 4*h*h/(1-h)*mp.quad(lambda s: w(s)*((1+mp.exp(-s))**h-1), [0, mp.inf])
    return A, R

for p in [2,3,7,20]:
    A,R = split(p)
    print(f"p={p}: |(4*2^-h + R) - pi_p| = {mp.nstr(abs(A+R-pi_cosh(p)),2)}")

# hypergeometric remark: g(h) = 4h int_0^1 v^(h-1)(1+v)^(-1-h) dv = 4*2F1(1+h,h;1+h;-1) = 4*2^(-h)
g = lambda h: 4*h*mp.quad(lambda v: v**(h-1)*(1+v)**(-1-h), [0,1])
for h in [mp.mpf('0.3'), mp.mpf('0.77')]:
    print(f"g({h}) - 4*2^(-h) = {mp.nstr(abs(g(h)-4*mp.power(2,-h)),2)}")

# J(m,k) = int_0^inf s^k log(1+e^-s)^m ds : key evaluations
def J(m,k): return mp.quad(lambda s: s**k*mp.log(1+mp.e**(-s))**m, [0, mp.inf])
print("J(1,0)=pi^2/12:", mp.nstr(abs(J(1,0)-pi**2/12),2))
print("J(1,1)=3zeta(3)/4:", mp.nstr(abs(J(1,1)-3*z3/4),2))
print("J(2,0)=zeta(3)/4:", mp.nstr(abs(J(2,0)-z3/4),2))
