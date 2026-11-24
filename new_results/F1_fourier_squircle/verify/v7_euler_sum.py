"""VERIFIER V7: Euler sums quoted in result.md section 4: sum H_k/k^2 = 2 zeta3, sum O_{2k}/k^2 = (7/4) zeta3,
O_{2k} = sum_{i<=k} 1/(2i-1).  Via H_k = int_0^1 (1-x^k)/(1-x) dx and O_{2k} = int_0^1 (1-x^{2k})/(1-x^2) dx:
   sum H_k/k^2 = int_0^1 (zeta2 - Li2(x))/(1-x) dx,   sum O_2k/k^2 = int_0^1 (zeta2 - Li2(x^2))/(1-x^2) dx."""
import mpmath as mp
mp.mp.dps = 40
z2 = mp.zeta(2)
s1 = mp.quad(lambda x: (z2 - mp.polylog(2, x)) / (1 - x), [0, mp.mpf(1)/2, 1])
s2 = mp.quad(lambda x: (z2 - mp.polylog(2, x * x)) / (1 - x * x), [0, mp.mpf(1)/2, 1])
print("sum H_k/k^2 - 2 zeta3 =", mp.nstr(s1 - 2 * mp.zeta(3), 3), ";  sum O_2k/k^2 - (7/4) zeta3 =", mp.nstr(s2 - mp.mpf(7)/4 * mp.zeta(3), 3))
