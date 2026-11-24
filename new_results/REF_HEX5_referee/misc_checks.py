"""Small independent checks quoted in the report."""
import mpmath as mp, sympy as sp
import direct5 as D, formula5 as FM, hpl
import exact_local_edge as XL
from local_examples import Zprinted

# 1. regular octagon P^(5)(1) (paper: 164.12775911330597065104567331245...)
mp.mp.dps = 60
m = 4
n = [(mp.cos(mp.pi * i / m), mp.sin(mp.pi * i / m)) for i in range(1, m + 1)]
R = 1 / (1 + mp.sqrt(2))
V = [(R * mp.cos(mp.pi / 2 + q * mp.pi / m), R * mp.sin(mp.pi / 2 + q * mp.pi / m)) for q in range(m)]
H = D.setup(V, [1] * m)
tot, per = D.all_derivs(H)
print('octagon: P(1) = %s (8(sqrt2-1)sqrt(2-sqrt2) = %s)' % (mp.nstr(tot[0], 30), mp.nstr(8 * (mp.sqrt(2) - 1) * mp.sqrt(2 - mp.sqrt(2)), 30)))
print('octagon: P5 = %s ; minus paper 164.12775911330597065104567331245 = %s' % (mp.nstr(tot[5], 50), mp.nstr(tot[5] - mp.mpf('164.12775911330597065104567331245'), 3)))
print('octagon: P5 minus OCT5D 120-digit value: %s' % mp.nstr(tot[5] - mp.mpf('164.12775911330597065104567331245436130507711746936109100978424924991432038932410908134739412239311309047813230825193813034902070966229339063253822702827777798757536177036'), 3))

# 2. local law on the paper's r = 2 edge: Z5 = -2079/250 sqrt5, Z23 = -477/250 sqrt5
V2, w2 = [(1, 0), (2, 2), (0, 1)], [1, 3, 1]
eds = FM.edges_data(V2, w2)
for q, ed in enumerate(eds):
    z5, z23 = Zprinted(ed)
    print('paper hexagon edge %d (r=%s): Z5 = %s, Z23 = %s' % (q, ed['r'], sp.nsimplify(z5), sp.nsimplify(z23)))

# 3. diamond: P5 = sqrt2(-90 z5 - 210 z2z3 + 396 z2^2 + 720 z3 - 480 z2)
E = XL.edge_exact(1, 1, 0, 2, -1, 1)          # edge (1,0)->(0,1): D=(-1,1), L2=2, q0=-1, V02=1
P5 = sp.expand(4 * E[5].subs({XL.m0: 0, XL.m1: 0, XL.mB: 0}))
z2 = sp.pi ** 2 / 6
target = sp.sqrt(2) * (-90 * sp.zeta(5) - 210 * z2 * sp.zeta(3) + 396 * z2 ** 2 + 720 * sp.zeta(3) - 480 * z2)
print('diamond exact P5 - printed =', sp.simplify(P5 - target))

# 4. the 13 HPLs at x = -1 against the exact reductions quoted in K1g step_51 / RUL1
mp.mp.dps = 60
a = hpl.atoms(-1)
l2 = mp.log(2); pi = mp.pi; z3 = mp.zeta(3); z5 = mp.zeta(5); L4 = mp.polylog(4, 0.5); L5 = mp.polylog(5, 0.5)
Z23 = lambda c: c * pi ** 2 * z3 / 6
ref = {
 'H1': -l2, 'Li2': -pi ** 2 / 12, 'Li3': -3 * z3 / 4, 'H011': z3 / 8, 'Li4': -7 * pi ** 4 / 720,
 'H0011': -pi ** 4 / 48 + mp.mpf(7) / 4 * z3 * l2 - pi ** 2 * l2 ** 2 / 12 + l2 ** 4 / 12 + 2 * L4,
 'H0111': -pi ** 4 / 90 + mp.mpf(7) / 8 * z3 * l2 - pi ** 2 * l2 ** 2 / 24 + l2 ** 4 / 24 + L4,
 'Li5': -mp.mpf(15) / 16 * z5,
 'H00011': -mp.mpf(29) / 32 * z5 + Z23(mp.mpf(1) / 2),
 'H00101': mp.mpf(41) / 32 * z5 - Z23(mp.mpf(5) / 8),
 'H00111': mp.mpf(33) / 32 * z5 + Z23(mp.mpf(1) / 2) - mp.mpf(7) / 8 * z3 * l2 ** 2 + pi ** 2 * l2 ** 3 / 18 - l2 ** 5 / 15 - 2 * L5 - 2 * L4 * l2,
 'H01011': -mp.mpf(73) / 64 * z5 - Z23(mp.mpf(3) / 2) + mp.mpf(7) / 4 * z3 * l2 ** 2 - pi ** 2 * l2 ** 3 / 9 + 2 * l2 ** 5 / 15 + 4 * L5 + 4 * L4 * l2,
 'H01111': z5 - mp.mpf(7) / 16 * z3 * l2 ** 2 + pi ** 2 * l2 ** 3 / 36 - l2 ** 5 / 30 - L5 - L4 * l2,
}
print('HPL(-1) vs quoted exact reductions, worst |diff| =', mp.nstr(max(abs(a[k] - ref[k]) for k in ref), 3))
