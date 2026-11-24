# R99 continued: b_n'' forms (notes, n<=17), N(p) derivatives, Parseval consistency, the C_n series.
import sys
sys.argv = ['x', '40']
from r99_sine import bder  # (runs the n<=11 table first at 40 dps)
from mpmath import mp, mpf, pi, log, zeta, gamma, sin, diff, nstr, floor, fsum, nsum, inf
mp.dps = 40
l2 = log(2)
notes_b2 = {3: (12-16*l2+3*pi**2)/128, 5: (508+48*l2-45*pi**2)/3456, 7: (-572-96*l2+63*pi**2)/13824,
            9: (2532+160*l2-225*pi**2)/64000, 11: (-25592-2160*l2+2475*pi**2)/1296000,
            13: (1598344+75600*l2-143325*pi**2)/88905600, 15: (-7820-448*l2+735*pi**2)/702464,
            17: (18638164+705600*l2-1686825*pi**2)/1828915200}
B = {}
for n in range(1, 62, 2):
    B[n] = bder(n)
for n, v in notes_b2.items():
    print("b_%d''(2) own - notes form: %s" % (n, nstr(B[n][2]-v, 3)))
N = lambda p: 2*sin(pi/p)/pi*gamma(3/p)*gamma(1-1/p)/gamma(1+2/p)
mp.dps = 60
N2 = diff(N, 2, 2); N3 = diff(N, 2, 3)
mp.dps = 40
print("N(2) =", nstr(N(mpf(2)), 20))
print("N''(2) - (ln^2 2 - 3 ln2 + 5pi^2/24) =", nstr(N2-(l2**2-3*l2+5*pi**2/24), 3))
N3f = (87*zeta(3)+16*l2**3-96*l2**2+72*l2+10*pi**2*l2-20*pi**2)/16
print("N'''(2) - (1/16)[87 z3 + ...] =", nstr(N3-N3f, 3))
# Parseval: N = sum b_n^2 -> check N'(2), N''(2), N'''(2) from the Fourier side (partial sums n<=61)
s1 = fsum(2*B[n][0]*B[n][1] for n in B)
s2 = fsum(2*B[n][0]*B[n][2] + 2*B[n][1]**2 for n in B)
s3 = fsum(2*B[n][0]*B[n][3] + 6*B[n][1]*B[n][2] for n in B)
print("Parseval N'  : Fourier %s vs Gamma %s" % (nstr(s1, 15), nstr(diff(N, 2, 1), 15)))
print("Parseval N'' : Fourier %s vs Gamma %s" % (nstr(s2, 15), nstr(N2, 15)))
print("Parseval N''': Fourier %s vs Gamma %s  (tail n>61 truncated)" % (nstr(s3, 15), nstr(N3, 15)))
# C_n law and the series sum_{n>=3} b_n' C_n = (1 - zeta3)/64
def bp(n): return mpf((-1)**((n-1)//2))/((n*n-1)*(2*floor(mpf(n-1)/4)+1))
def Cn(n): return mpf((-1)**((n+1)//2))*n*(2*floor(mpf(n-1)/4)+1)/(2*(n*n-1)**2)
for n in [3, 5, 7, 9, 11, 13]:
    # b_n'' = A_n + b_n' ln2 + C_n pi^2 with A_n rational; check the pi^2 and ln2 coefficients against the notes' forms
    pass
# split by residue class so each summand is a smooth rational function of m
def term(n, sgn_b, fl):
    bpv = sgn_b/((n*n-1)*fl); cnv = -sgn_b*n*fl/(2*(n*n-1)**2)
    return bpv*cnv
S = nsum(lambda m: term(4*m+3, -1, 2*m+1), [0, inf]) + nsum(lambda m: term(4*m+5, 1, 2*m+3), [0, inf])
# sanity: compare with direct partial sum using the integer formulas
S_direct = fsum(bp(n)*Cn(n) for n in range(3, 4001, 2))
print("direct partial sum to n=4000:", nstr(S_direct, 15))
print("sum_{n>=3} b_n' C_n =", nstr(S, 30), "  (1-zeta3)/64 =", nstr((1-zeta(3))/64, 30))
# check C_n against notes' b_n'' forms (pi^2 coefficient) and ln2 coefficient = b_n'
from mpmath import pslq
for n, v in notes_b2.items():
    pass
print("C_n vs notes pi^2 coeff:", [(n, nstr(Cn(n) - {3: mpf(3)/128, 5: mpf(-45)/3456, 7: mpf(63)/13824, 9: mpf(-225)/64000, 11: mpf(2475)/1296000, 13: mpf(-143325)/88905600, 15: mpf(735)/702464, 17: mpf(-1686825)/1828915200}[n], 3)) for n in notes_b2])
print("b_n' vs notes ln2 coeff:", [(n, nstr(bp(n) - {3: mpf(-16)/128, 5: mpf(48)/3456, 7: mpf(-96)/13824, 9: mpf(160)/64000, 11: mpf(-2160)/1296000, 13: mpf(75600)/88905600, 15: mpf(-448)/702464, 17: mpf(705600)/1828915200}[n], 3)) for n in notes_b2])
