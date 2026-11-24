import mpmath as mp, json
mp.mp.dps = 100
pi4 = mp.pi**4
fit = [mp.mpf(s) for s in json.load(open('fit_w2_w12.json'))]
# exact residuals: w2 = z2/4, w3 = z3/4, w4 = -19 pi^4/2880, w5 = -z2 z3/2 - 3 z5/16
tests = [
    ("w2 - z2/4", fit[0] - mp.zeta(2)/4),
    ("w3 - z3/4", fit[1] - mp.zeta(3)/4),
    ("w4 - (-19 pi^4/2880)", fit[2] + 19*pi4/2880),
    ("w5 - (-z2z3/2 - 3z5/16)", fit[3] + mp.zeta(2)*mp.zeta(3)/2 + 3*mp.zeta(5)/16),
    ("w4/gamma guard: (w4+19pi^4/2880)/gamma", (fit[2] + 19*pi4/2880)/mp.euler),
]
for name, r in tests:
    print(f"{name} = {mp.nstr(r, 4)}")
