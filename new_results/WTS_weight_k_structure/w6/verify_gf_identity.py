"""verify step: exact check of the generating-function identity quoted in identify_w6.md:
   sum_m eps_m A_{m+1} z^m/m!  =  (1 - sinh z)/cosh^2 z,  eps_m = +,-,-,+ repeating, A_n Euler zigzag numbers,
   and of the rational prefactors eps_m A_{m+1}/(2 m!) for m = 0..6."""
import sympy as sp
z = sp.symbols('z')
N = 14
ser = sp.series((1 - sp.sinh(z)) / sp.cosh(z) ** 2, z, 0, N).removeO()
A = [sp.Integer(1)] + [sp.Integer(x) for x in (1, 1, 2, 5, 16, 61, 272, 1385, 7936, 50521, 353792, 2702765, 22368256, 199360981)]
eps = lambda m: [1, -1, -1, 1][m % 4]
# independent zigzag via tan+sec series: sum A_n z^n/n!
zz = sp.series(sp.sec(z) + sp.tan(z), z, 0, N + 1).removeO()
Aser = [zz.coeff(z, n) * sp.factorial(n) for n in range(N + 1)]
print("zigzag numbers from sec+tan series match the table used:", [Aser[n] for n in range(N + 1)] == A[:N + 1])
ok = all(sp.simplify(ser.coeff(z, m) - eps(m) * Aser[m + 1] / sp.factorial(m)) == 0 for m in range(N - 1))
print("series of (1-sinh z)/cosh^2 z equals sum eps_m A_(m+1) z^m/m! through z^%d:" % (N - 2), ok)
print("prefactors eps_m A_(m+1)/(2 m!), m=0..6:", [sp.Rational(eps(m) * Aser[m + 1], 2 * sp.factorial(m)) for m in range(7)])
