# P2N4R referee: recompute every constant of the proof of Lemma g3:prop:n4 in exact rational arithmetic.
# Only inequalities between rationals and exact sympy reals (log 9, e^{1/4}, ...), all decided by sympy exactly.
import sympy as sp
from sympy import Rational as Q
out = []
def p(*s):
    print(*s); out.append(' '.join(map(str, s)))
def le(a, b):           # decide a <= b exactly (sympy, high-precision interval evaluation of a-b)
    d = sp.N(sp.sympify(b) - sp.sympify(a), 60); return bool(d >= 0)
beta = Q(9, 8)
L0 = Q(21972, 10000)
p("log 9 >= 2.1972:", le(L0, sp.log(9)), "; log 9 =", sp.N(sp.log(9), 12))
# --- chain for sup|H| <= 3.1 on |eps|<=1/4 (Lemmas 2.2, 2.6) ---
sumtk = sp.exp(Q(1, 2))/4*sp.zeta(Q(3, 2))
p("sum|t_k| <= e^{1/2}/4 zeta(3/2) =", sp.N(sumtk, 10), "-> |F| <=", sp.N(1+sumtk, 10))
pis = sp.sinh(sp.pi/4)/(sp.pi/4) - 1
p("|pi e/sin pi e| <= 1/(1-(sinh(pi/4)/(pi/4)-1)) =", sp.N(1/(1-pis), 10))
gsup = (1+sumtk)/(1-pis)
p("|g| <=", sp.N(gsup, 10), "<= 2.4:", le(gsup, Q(24, 10)))
Hsup = gsup*sp.exp(Q(1, 4))
p("|H| <= |g| e^{|a1|/4} <=", sp.N(Hsup, 10), "<= 3.1:", le(Hsup, Q(31, 10)), "; with printed 2.4:", sp.N(Q(24,10)*sp.exp(Q(1,4)), 8))
# --- Step 0 ---
p("beta^5 =", float(beta**5), "< 1.8021:", beta**5 < Q(18021, 10000), "; beta^6 =", float(beta**6), "< 2.0273:", beta**6 < Q(20273, 10000))
s1 = (sp.log(9)+1)/9; s2 = (sp.log(9)+1)**2/9
p("(log9+1)/9 =", sp.N(s1, 8), "< 0.3553:", le(s1, Q(3553, 10000)), "; (log9+1)^2/9 =", sp.N(s2, 8), "< 1.1359:", le(s2, Q(11359, 10000)))
# monotonicity claims of (iii): d/dn[(log n+1)/n] = -log n/n^2 <0 ; d/dn[(log n+1)^2/n] = (log n+1)(1-log n)/n^2 <0 for n>e
n = sp.symbols('n', positive=True)
p("d/dn (log n+1)/n =", sp.simplify(sp.diff((sp.log(n)+1)/n, n)), "; d/dn (log n+1)^2/n =", sp.factor(sp.simplify(sp.diff((sp.log(n)+1)**2/n, n))))
p("disc |e-e0|<=1/10, e0<=1/8 inside |e|<=0.225<1/4:", Q(1,8)+Q(1,10) == Q(225,1000))
A = Q(215, 100); a3m = Q(31, 10)*5**3; b4m = Q(31, 10)*5**4; K5 = Q(31, 10)*10**5
p("|a3| <= 3.1*5^3 =", a3m, "; |b4| <= 3.1*5^4 =", b4m, "; |H^(5)|/5! <= 3.1*10^5 =", K5)
p("A = zeta(2)+1/2 < 2.15:", le(sp.zeta(2)+Q(1, 2), A))
# --- Step 5 constants ---
m2 = Q(9, 8)/(1-Q(3, 2)*Q(1, 9)); m3 = Q(4, 3)/(3-Q(4, 9)); m4 = Q(9, 31)*beta
p("m=2:", m2, "=1.35:", m2 == Q(135, 100), "; m=3:", m3, "=12/23:", m3 == Q(12, 23), "; m>=4:", m4, "=81/248:", m4 == Q(81, 248))
p("(8m-1)/9 at m=4 =", Q(31, 9), "; check m-(m+1)/9 = (8m-1)/9")
# --- Step 7: C0, C1, C2 with the printed rounded numbers 0.3553, 1.1359 ---
S1, S2 = Q(3553, 10000), Q(11359, 10000)
b5, b6 = beta**5, beta**6
pieces = {
 'e1 (r2,r3,r4)': (2*A*b5*4, 2*A*b5*6, 2*A*b5*2),
 'e1 (eps^5 lam^3 term)': (0, 0, A/3*(5*S1*b6 + 3 + S2*b6)),
 'e2': (2*a3m*6*b5, 2*a3m*4*b5, 2*a3m*b5/2),
 'e3': (2*b4m*4*b5, 2*b4m*b5, 0),
 'P4': (2*K5*b5, 0, 0),
 'e5': (Q(144, 7), 0, 0),
 'e6': (2*(m2+m3+m4), 0, 0),
}
C0 = sum(v[0] for v in pieces.values()); C1 = sum(v[1] for v in pieces.values()); C2 = sum(v[2] for v in pieces.values())
for k, v in pieces.items(): p("  %-24s" % k, [float(t) for t in v])
p("C0 =", float(C0), "<= 1.154e6:", C0 <= Q(1154000)); p("C1 =", float(C1), "<= 1.262e4:", C1 <= 12620); p("C2 =", float(C2), "<= 721:", C2 <= 721)
p("P4 share of C0: 2*3.1e5*beta^5 =", float(2*K5*b5), "< 1.118e6:", 2*K5*b5 < Q(1118000))
kappa = 1 + 1/L0
p("kappa = 1+1/2.1972 =", float(kappa), "< 1.4552:", kappa < Q(14552, 10000))
Cp = C0/L0**2 + C1*kappa/L0 + C2*kappa**2
p("C' = C0/L0^2 + C1 kappa/L0 + C2 kappa^2 =", float(Cp), "(exact rational, L0=2.1972) ; < 2.49e5:", Cp < 249000, "; <= 2.5e5:", Cp <= 250000)
C4 = A/3 + Cp/L0
p("C4 = 2.15/3 + C'/L0 =", float(C4), "; < 1.14e5:", C4 < 114000, "; <= 1.2e5:", C4 <= 120000)
# the same with the stated rounded C0,C1,C2 (as the LaTeX text does)
Cp2 = Q(1154000)/L0**2 + Q(12620)*Q(14552, 10000)/L0 + 721*Q(14552, 10000)**2
p("C' from the printed roundings (1.154e6, 1.262e4, 721, kappa 1.4552) =", float(Cp2), "< 2.49e5:", Cp2 < 249000)
p("C4 from that:", float(A/3 + Cp2/L0))
# slack: what the bound would be at the true sup|H| (not needed)
open('ref_constants.log', 'w').write('\n'.join(out) + '\n')
