# Planner check of COR7's c = 1/3 closed forms, transcribed from result.md text (not from its JSON).
import sys, mpmath as mp
D = int(sys.argv[1]) if len(sys.argv) > 1 else 400
mp.mp.dps = D + 30
d = "<project>/research/08_theory/COR7_classical_corners/"
rd = lambda f: mp.mpf(open(d + f).read().strip().split()[0])
M3q = rd("c13_M3_quad_1100.txt"); W4h = rd("c13_W4_hpl_1150.txt"); W3h = rd("c13_W3_hpl_1150.txt"); W2h = rd("c13_W2_hpl_1150.txt")
L2, L3, pi, z3 = mp.log(2), mp.log(3), mp.pi, mp.zeta(3)
Li = mp.polylog
l13, l12m = Li(2, mp.mpf(1)/3), None
W2 = -Li(2, mp.mpf(1)/3) - Li(2, -mp.mpf(1)/3) + pi**2/12
W3 = (-26*z3 - 87*Li(3, mp.mpf(1)/3) - 130*Li(3, -mp.mpf(1)/2) + mp.mpf(21)/4*pi**2*L2 - mp.mpf(43)/12*pi**2*L3
      + 67*l13*L2 - 45*l13*L3 + mp.mpf(65)/3*L2**3 - 65*L2**2*L3 + mp.mpf(109)/2*L2*L3**2 - mp.mpf(37)/3*L3**3)
t3 = Li(3, mp.mpf(1)/3); tm = Li(3, -mp.mpf(1)/2)
W4 = (-mp.mpf(22)/3*Li(4, -2) - mp.mpf(17)/2*Li(4, mp.mpf(1)/3) + mp.mpf(20)/3*Li(4, mp.mpf(2)/3) + mp.mpf(110)/9*Li(4, mp.mpf(1)/2)
      + mp.mpf(25)/48*Li(4, mp.mpf(1)/9)
      - 12*t3*L2 + mp.mpf(28)/3*t3*L3 - mp.mpf(16)/3*tm*L2 + mp.mpf(16)/3*tm*L3
      + mp.mpf(14)/3*z3*L2 - mp.mpf(25)/9*z3*L3
      - 4*l13**2 + mp.mpf(13)/3*l13*L2**2 - 6*l13*L2*L3 + mp.mpf(2)/3*l13*L3**2 + mp.mpf(11)/18*pi**2*l13
      - mp.mpf(31)/27*pi**2*L2**2 + mp.mpf(31)/72*pi**2*L3**2 - mp.mpf(271)/1296*pi**4
      + mp.mpf(46)/27*L2**4 - mp.mpf(62)/9*L2**3*L3 + mp.mpf(61)/6*L2**2*L3**2 - mp.mpf(17)/3*L2*L3**3 + mp.mpf(103)/144*L3**4)
for nm, a, b in [("W2", W2, W2h), ("W3", W3, W3h), ("W4", W4, W4h), ("M3 total vs quadrature", W2 + W3 + W4, M3q)]:
    print(f"{nm}: |closed form - value| = {mp.nstr(abs(a - b), 5)}")
print("dps", D)
