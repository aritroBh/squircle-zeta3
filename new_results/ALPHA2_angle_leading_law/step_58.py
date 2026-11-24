from mpmath import mp, mpf, quad, nstr, pi, sin, cos, legenp
mp.dps = 60
def Jfull_careful(s, alpha):
    # split the interval finely near endpoints
    pts = [0]
    x = mpf('1e-3')
    while x < (pi-alpha)/2:
        pts.append(x); x *= 8
    mid = pts + [pi-alpha-p for p in reversed(pts[:-1])] + [pi-alpha]
    mid = sorted(set(mid))
    return quad(lambda u: sin(u)**(s-1)*sin(alpha+u)**(-s), mid)
for alpha in [pi/2, 2*pi/3, pi/6]:
    ca = cos(alpha)
    for s in [mpf('0.9'), mpf('0.95'), mpf('0.1')]:
        a = Jfull_careful(s, alpha); b = pi*legenp(s-1, 0, -ca, type=2)/sin(pi*s)
        print(f"alpha={nstr(alpha,6)}, s={nstr(s,4)}: |diff|={nstr(abs(a-b),3)}")
