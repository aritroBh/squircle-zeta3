import mpmath as mp, time
import cx6r_common as C
print("CX6R step_16: sec.3 generating-function identity + (15)/(16) branch checks")
mp.mp.dps = 80
def M_orig(eps):
    def f(t):
        a = t**(eps-1); b = (1-t)**(eps-1)
        return a + b - a*mp.sqrt(1+(b/a)**2)
    return -C.glint(f, [[mp.mpf('1e-30'), mp.mpf('0.25')],[mp.mpf('0.25'), mp.mpf('0.5')]], 6)
def Phi(eps, s):
    return C.glint(lambda t: t**(s*(1-eps))*(1+t)**(-1-eps), [[0,0.5],[0.5,1]], 5)
def J(eps): return (1-mp.mpf(2)**(-eps))/eps if eps != 0 else mp.log(2)
print("(a) M(eps) = -J(eps) + (2/pi) Im int_Gamma H K Phi(eps,s) ds")
for eps in [mp.mpf('0.05'), mp.mpf('-0.1'), mp.mpf('0.2')]:
    I = C.glint(lambda r: C.H(mp.mpf('0.5')+r*C.w4)*C.K(mp.mpf('0.5')+r*C.w4)*Phi(eps, mp.mpf('0.5')+r*C.w4)*C.w4, C.GAMMA_SEGS, 5)
    val = -J(eps) + 2/mp.pi*mp.im(I)
    print(f"  eps={mp.nstr(eps,5)}: |M_orig - (-J + 2/pi Im int)| = {mp.nstr(abs(val-M_orig(eps)),5)}")
print("(b) i<=k/2 saddle comparisons at k=99 (tightest first):")
k = 99; N = 100
rows = []
for i in range(1, k//2+1):
    lhs = mp.sqrt((3*mp.pi/4)*(i+1)*(k-i))
    rows.append((lhs - 2*mp.sqrt(N), i))   # must be >= 0
rows.sort()
print("   min sqrt((3pi/4)(i+1)(k-i)) - 2 sqrt(N) over 1<=i<=49:", mp.nstr(rows[0][0],10), "at i =", rows[0][1])
print("   i=0: sqrt((3pi/2)k) - 2 sqrt(N) =", mp.nstr(mp.sqrt(3*mp.pi/2*k)-2*mp.sqrt(N),10))
print("(c) i>k/2 branch (16) at k=99:")
x0 = 2*mp.exp(mp.pi/4)
s_hi = sum(x0**i/mp.factorial(i) for i in range(k//2+1, k+1))
bound = mp.exp(-k/mp.mpf(2))*mp.exp(mp.e*x0)
print(f"   sum_{{i>k/2}} x^i/i! = {mp.nstr(s_hi,10)}   <= e^-k/2 e^{{ex}} = {mp.nstr(bound,10)}   ratio {mp.nstr(s_hi/bound,6)}")
print("   e^{x} =", mp.nstr(mp.exp(x0), 15), "  e^{e x} =", mp.nstr(mp.exp(mp.e*x0), 15))
print("(d) (17) exponent: sqrt((9pi/4)k) - 2 sqrt(N) =", mp.nstr(mp.sqrt(9*mp.pi/4*k)-2*mp.sqrt(N),10))
