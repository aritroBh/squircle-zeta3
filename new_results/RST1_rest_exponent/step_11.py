"""RST1 step_11: large-k behavior of bound/target; piece decomposition."""
import mpmath as mp
from rst1_common import *
print("bound(k)/target(k) for large k:")
for k in [99,200,400,1000,3000,10000,30000,100000]:
    print(f"  k={k:7d}: {mp.nstr(rest_bound(k)/target(k),8)}")
print()
print("piece decomposition (normalized to target) at k=99,400,10000:")
for k in [99,400,10000]:
    i0 = FH_gauss(k,2); i1=2*FH_gauss(k-1,2)
    S = mpf(0)
    for i in range(2, k//2+1): S += 2**i/mp.factorial(i)*FH_box(k-i,i+1)
    T = mpf(0)
    for i in range(k//2+1, k+1): T += 2**i/mp.factorial(i)*FH_box(k-i,i+1)
    t=target(k)
    print(f" k={k}: q2={(8/pi)*(i0+i1)/t!s:>12.6}  mid={(8/pi)*S/t!s:>12.6}  tail={(8/pi)*T/t!s:>12.6}  T2={T2_bound(k)/t!s:>12.6}  T3J={T3J_bound(k)/t!s:>12.6}")
