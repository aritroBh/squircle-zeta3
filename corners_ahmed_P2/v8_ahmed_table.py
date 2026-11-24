"""V8. Reproduce Ahmed's Table IV (empirical constants hat C_n(theta) at n = 10^6) from the exact formula,
and compare with the second-order prediction hat C_n = C - [z (log n - a1)^2 + Q2]/(n-1) + O(log^3 n/n^2)."""
import mpmath as mp
from common import *
from closed_forms import closed_a
setdps(40)
ahmed = {'pi/6': '0.0045781', 'pi/4': '0.0225492', 'pi/3': '0.0684278', 'pi/2': '0.3067523',
         '2pi/3': '0.8066971', '3pi/4': '1.1442509', '5pi/6': '1.5036576'}
angs = {'pi/6': mp.pi/6, 'pi/4': mp.pi/4, 'pi/3': mp.pi/3, 'pi/2': mp.pi/2, '2pi/3': 2*mp.pi/3,
        '3pi/4': 3*mp.pi/4, '5pi/6': 5*mp.pi/6}
n = mp.mpf(10)**6
print(f"{'theta':>6} {'C(theta)':>12} {'hatC_n (exact formula)':>26} {'Ahmed Tab.IV':>13} {'2nd-order pred.':>18}")
for k, th in angs.items():
    z = mp.sin(th/2)**2
    A = closed_a(th)
    C = 2*z + 2*(1 - z)*mp.log(1 - z)
    hat = (n - 1)*(2 - L_exact(n, th)) - 2*z*mp.log(n)
    Q2 = 2*z*(A[1] - (1 - z))
    pred = C - (z*(mp.log(n) - A[0])**2 + Q2)/(n - 1)
    print(f"{k:>6} {mp.nstr(C,8):>12} {mp.nstr(hat,15):>26} {ahmed[k]:>13} {mp.nstr(pred,15):>18}   |hat-pred|={mp.nstr(abs(hat-pred),2)}")
