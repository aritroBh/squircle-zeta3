# Closed forms of thm:S at full precision (s4_S.py built log2, log(1+sqrt2), sqrt2-1, arctan(1/sqrt2)
# at mpmath's default 15 digits, so its closed-form comparisons are only good to ~1e-15).
# Compare against the numerical values printed in s4_S.out.
import re
from mpmath import mp, mpf, sqrt, log, pi, atan, polylog, zeta, nstr
mp.dps = 50
txt = open('s4_S.out').read()
num = lambda pat: mpf(re.search(pat, txt).group(1))
Q0n = num(r"Q0: integral (\S+)"); Q1n = num(r"Q1 \(derivative under the integral\) (\S+)")
M0n = num(r"== M0 (\S+)"); M1n = num(r"== M1 (\S+)"); M2n = num(r"== M2 (\S+)")
s1n = num(r"== s1 (\S+)"); s2n = num(r"== s2 (\S+)"); s3n = num(r"== s3 (\S+)")
lam = log(2); L = log(1 + sqrt(2)); b = sqrt(2) - 1; ac = atan(1 / sqrt(2))
Q0c = (-pi**2/8 + sqrt(2)/16*pi**2 + 3*sqrt(2)/2*pi*ac - 3*lam**2 - 9*sqrt(2)*L*log(3)
       + (6 + 6*sqrt(2))*L**2 - 6*sqrt(2)*polylog(2, b) + 9*sqrt(2)/4*polylog(2, -b**4))
M0c = lam - sqrt(2)*L
M1c = (19*sqrt(2)/48 - mpf(1)/24)*pi**2 + lam**2/2 + 2*sqrt(2)*lam*L - (1 + 5*sqrt(2)/2)*L**2 - 5*sqrt(2)*polylog(2, b)
s2c = (-5*pi**2 + 37*sqrt(2)/2*pi**2 - 12*sqrt(2)*pi*ac + 48*lam**2 - 96*L**2 - 168*sqrt(2)*L**2
       + 96*sqrt(2)*lam*L + 72*sqrt(2)*log(3)*L - 192*sqrt(2)*polylog(2, b) - 18*sqrt(2)*polylog(2, -b**4))
out = []
out.append("Q0 closed (paper, thm:S) = %s ; my integral - closed = %s" % (nstr(Q0c, 40), nstr(Q0n - Q0c, 3)))
out.append("M0 closed - numeric = %s" % nstr(M0c - M0n, 3))
out.append("M1 closed (thm:M1) - numeric = %s" % nstr(M1c - M1n, 3))
out.append("s1: 48*M0closed = %s ; numeric s1 diff %s" % (nstr(48*M0c, 30), nstr(48*M0c - s1n, 3)))
s2fromclosed = 48*M1c - 4*pi**2 - 8*Q0c
out.append("s2 from closed M1,Q0 = %s ; T6 closed form - that = %s ; numeric s2 - T6 = %s" % (nstr(s2fromclosed, 30), nstr(s2c - s2fromclosed, 3), nstr(s2n - s2c, 3)))
out.append("s3 numeric = %s (printed 48.2874324860...)" % nstr(s3n, 30))
out.append("Q1 numeric = %s (printed 6.5043255737...)" % nstr(Q1n, 30))
out.append("printed decimals: s1 %s s2 %s" % (nstr(48*M0c, 12), nstr(s2c, 12)))
print("\n".join(out)); open('s4b_constants.out', 'w').write("\n".join(out) + "\n")
