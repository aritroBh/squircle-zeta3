"""Certified (arb) constant C_U in |U_k| <= C_U * exp(-(3/2) sqrt(pi (N-3))), all N >= 100.
U_k = rest - S_k = -J_k + (2/pi) Im{ I1 + I2 + I3 - int_0^{1/2} H E Q_k - int_0^{1/2} H E A2_k }.
Uses CX6 (13),(14),(16),(17) and the new bound (B3) on A_k^{>=3}. Each piece is written as
(const) * exp(-c sqrt(N - n0)) with c sqrt(N-n0) >= (3/2) sqrt(pi (N-3)), proved by hand in result.md."""
from flint import arb, ctx
ctx.prec = 300
pi = arb.pi(); a = pi/(2*arb(2).sqrt()); r2 = arb(2).sqrt()
def Fpref(q, lam):  # F(j,q;lam) <= Fpref * exp(-sqrt(3 lam j q / sqrt2))
    return 4/lam*(lam/2 + lam*q/r2).exp()
x = 2*(pi/4).exp()
# I1: (2/pi)|I1| <= (2/pi)*8*F(k-1,2;3a); exponent sqrt(9 pi (N-2)/2) >= (3/2)sqrt(pi(N-3))
c1 = 2/pi*8*Fpref(2, 3*a)
# I2: (2/pi)|I2| <= (2/pi)*4*[2F(k,3;a)+4F(k-1,3;a)+2 sum_{i>=2} 2^i/i! F(k-i,i+1;a)]
#   F(k,3), F(k-1,3): exp(-sqrt(9pi(N-2)/4)) ; i in [2,k/2]: (k-i)(i+1)>=3(k-2): exp(-sqrt(9 pi (N-3)/4));
#   prefactor of F(k-i,i+1;a) = (4/a) e^{a/2} e^{pi(i+1)/4};  sum_i 2^i/i! e^{pi i/4} <= e^{x}
#   i > k/2 : drop the saddle factor (it is <=1), sum_{i>k/2} x^i/i! <= e^{e x} e^{-k/2}  (CX6 (16))
pa = 4/a*(a/2).exp()*(pi/4).exp()
c2 = 2/pi*4*(2*Fpref(3, a) + 4*Fpref(3, a) + 2*pa*x.exp() + 2*pa*(arb(1).exp()*x).exp())
# I3: CX6 (17): (2/pi)|I3| <= T0 exp(-sqrt(9 pi k/4)),  T0 = 32/(3 pi a) e^{3a/2+3pi/4}
T0 = 32/(3*pi*a)*(3*a/2 + 3*pi/4).exp()
# elementary: J_k <= (log2)^N/N! <= e^{-(N-1)/2};  segments <= (2/(pi k)) 3^{-k} + (2/pi)*3*2^{-2k}/(k-1) <= e^{-(N-1)/2}
c_el = arb(2)
CU = c1 + c2 + T0 + c_el
print("c1 (I1)   =", c1.str(15)); print("c2 (I2)   =", c2.str(15)); print("T0 (I3)   =", T0.str(15))
print("C_U       =", CU.str(15), " upper:", CU.upper())
# check e^{-k/2} <= e^{-(3/2)sqrt(pi(N-3))} and e^{-(N-1)/2} likewise at N=100 (monotone beyond, see result.md)
N = arb(100)
print("at N=100: k/2 =", ((N-1)/2).str(8), " (3/2)sqrt(pi(N-3)) =", (arb(3)/2*(pi*(N-3)).sqrt()).str(8))
print("sqrt(9pi/4) =", (9*pi/4).sqrt().str(12), " sqrt(2pi) =", (2*pi).sqrt().str(12), " sqrt(3pi) =", (3*pi).sqrt().str(12))
