# Planner check of COR9 (independent of COR9's code): weight split W2/W3/W4 at the three corners from G1's proved B3 JSON
# with the planner's own HPL series evaluator, then the closed forms transcribed from COR9 result.md.
import json, re, sys, mpmath as mp
from fractions import Fraction as Fr
D = int(sys.argv[1]) if len(sys.argv) > 1 else 120
mp.mp.dps = D + 20
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
C9 = "<project>/research/08_theory/COR9_golden_corner_closed_forms/"
def hpl(w, x, N):
    a = w[-1]; s = [mp.mpf(0)] * (N + 1)
    for m in range(1, N + 1): s[m] = mp.mpf(1) / m if a == 1 else mp.mpf((-1) ** (m + 1)) / m
    for a in reversed(w[:-1]):
        t = [mp.mpf(0)] * (N + 1)
        if a == 0:
            for m in range(1, N + 1): t[m] = s[m] / m
        else:
            acc = mp.mpf(0)
            for m in range(0, N): acc = acc * a + s[m]; t[m + 1] = acc / (m + 1)
        s = t
    return mp.fsum(s[m] * x ** m for m in range(1, N + 1))
def terms(f):
    out = []
    for nm, v in json.load(open(f)):
        q = mp.mpf(Fr(v).numerator) / Fr(v).denominator
        if nm == '1': out.append(('1', (), 0, q)); continue
        m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
        out.append((m.group(2), tuple(int(t) for t in m.group(3).split(',')), int(m.group(1)), q))
    return out
def split(c):
    N = int((D + 25) / -mp.log10(c)) + 30; cache = {}
    def tv(kind, w, e):
        if kind == '1': return mp.mpf(1)
        if (w, 1) not in cache: cache[(w, 1)] = hpl(w, c, N)
        f = c ** e * cache[(w, 1)]
        if kind == 'He':
            if (w, -1) not in cache: cache[(w, -1)] = hpl(w, -c, N)
            f += (-c) ** e * cache[(w, -1)]
        return f
    W = {}
    for fname, zw, zval in [("hpl_fit_k3_(0,0,0).json", 0, mp.mpf(1)), ("hpl_fit_k3_(1,0,0).json", 2, mp.zeta(2)), ("hpl_fit_k3_(0,1,0).json", 3, mp.zeta(3))]:
        for kind, w, e, q in terms(G1 + fname):
            W[len(w) + zw] = W.get(len(w) + zw, 0) + q * zval * tv(kind, w, e)
    return W
Li = lambda n, x: mp.polylog(n, x); pi2 = mp.pi**2; z3 = mp.zeta(3); lam = mp.log(2); Q = mp.mpf
s2 = mp.sqrt(2); b = s2 - 1; L2 = mp.log(1 + s2)
s5 = mp.sqrt(5); phi = (1 + s5)/2; L = mp.log(phi)
cf = {}
# ---- (sqrt2-1)^2
cf['b2sq'] = (b**2, {
 2: -Li(2, b**4)/2 + pi2/12,
 3: (-Q(273)/2 - Q(1655)/16*s2)*z3 + (-1166 - 898*s2)*Li(3, b) + (562 + 422*s2)*Li(3, -b) + (1544 + 1180*s2)*Li(3, 1/s2)
    + (-314 - 236*s2)*L2*Li(2, b) + (392 + 295*s2)*lam*Li(2, b) + (-Q(747)/8 - Q(287)/4*s2)*pi2*L2 + (Q(249)/4 + Q(767)/16*s2)*pi2*lam
    + (-Q(193)/6 - Q(295)/12*s2)*lam**3 + (-193 - Q(295)/2*s2)*lam**2*L2 + (272 + Q(413)/2*s2)*lam*L2**2 + (Q(56)/3 + Q(47)/3*s2)*L2**3,
 4: (-478 - 153*s2)*Li(4, b) + (-470 - 153*s2)*Li(4, -b) + (Q(83)/2 + Q(59)/4*s2)*Li(4, -b**2) + (4 - 4*s2)*Li(4, 1/s2)
    + (Q(103)/2 + Q(41)/4*s2)*Li(4, Q(1)/2) + (16 + 14*s2)*Li(4, 1 - 1/s2) + (116 + 30*s2)*Li(4, 2 - s2)
    + (12 + 8*s2)*L2*Li(3, b) + (-8 - 4*s2)*L2*Li(3, 1/s2) + (-Q(9)/4 - Q(3)/2*s2)*L2*z3
    + (-14 - 10*s2)*lam*Li(3, b) + (8 + 4*s2)*lam*Li(3, 1/s2) + (Q(17)/4 + Q(21)/8*s2)*lam*z3
    + (-15 - 10*s2)*Li(2, b)**2 + (-6 - 4*s2)*Li(2, b)*L2**2 + (-8 - 5*s2)*Li(2, b)*lam*L2 + (5 + 3*s2)*Li(2, b)*lam**2 + (Q(13)/6 + Q(17)/12*s2)*pi2*Li(2, b)
    + (-Q(361)/320 - Q(1681)/5760*s2)*mp.pi**4 + (Q(7)/3 + Q(5)/4*s2)*pi2*L2**2 + (Q(25)/6 + Q(11)/24*s2)*pi2*lam*L2 + (-Q(43)/12 - Q(3)/4*s2)*pi2*lam**2
    + (-Q(1)/4 - Q(25)/24*s2)*L2**4 + (-Q(52)/3 - Q(9)/2*s2)*L2**3*lam + (Q(61)/4 + Q(35)/8*s2)*L2**2*lam**2 + (-Q(9)/2 - Q(5)/8*s2)*L2*lam**3 + (Q(235)/96 + Q(13)/64*s2)*lam**4})
# ---- 1/phi
l23 = Li(2, phi**-3)
cf['gol1'] = (1/phi, {
 2: pi2/20 + L**2/2,
 3: (Q(4)/5 - s5/10)*z3 - Li(3, 1/phi) + (-Q(41)/240 - s5/16)*pi2*L + (-Q(7)/24 - 5*s5/24)*L**3,
 4: (-Q(5)/2 + Q(11)/6*s5)*Li(4, 1/(2*phi**2)) + (Q(113)/16 - Q(187)/96*s5)*Li(4, phi**-2) + (Q(3)/10 + Q(19)/30*s5)*Li(4, Q(1)/2)
    + (-Q(7)/2 + Q(3)/2*s5)*Li(4, 2/phi**2) + (-2 - Q(2)/3*s5)*Li(4, phi/2) + (-Q(29)/72 + s5/12)*Li(4, phi**-6) + (Q(35)/18 + Q(7)/6*s5)*Li(4, -phi**-3)
    + (Q(339)/40 + Q(113)/40*s5)*L*z3 + (-1 - 2*s5)*L*Li(3, 1/phi) + (-5 + Q(2)/3*s5)*L*Li(3, phi**-3) + (-5 - s5)*L*Li(3, 2/phi**2)
    + (-Q(7)/30 + s5/60)*pi2*l23 + l23**2 + (Q(25)/4 + Q(3)/4*s5)*l23*L**2 - 6*l23*L*lam
    + (Q(64)/2025 - s5/5400)*mp.pi**4 + (-Q(77)/40 - Q(67)/180*s5)*pi2*L**2 + (Q(1)/5 + Q(31)/180*s5)*pi2*L*lam + (Q(3)/10 - s5/30)*pi2*lam**2
    + (Q(411)/16 + Q(2041)/288*s5)*L**4 + (-24 - Q(143)/18*s5)*L**3*lam + (6 + Q(5)/3*s5)*L**2*lam**2 + (-Q(1)/2 + Q(13)/18*s5)*L*lam**3 + (-Q(7)/40 + Q(3)/40*s5)*lam**4})
# ---- sqrt5-2 = 1/phi^3
cf['gol3'] = (phi**-3, {
 2: -Li(2, phi**-6)/2 + pi2/12,
 3: (-Q(2029)/40 - Q(1061)/40*s5)*z3 + (-19 - Q(17)/2*s5)*Li(3, 1/phi) + (Q(35)/3 + 5*s5)*Li(3, phi**-3) + (62 + 32*s5)*Li(3, 2/phi**2)
    + (Q(151)/10 + Q(116)/15*s5)*pi2*L + (-Q(47)/6 - 4*s5)*pi2*lam + (-48 - 24*s5)*l23*L + (32 + 16*s5)*l23*lam
    + (-Q(3905)/12 - Q(2003)/12*s5)*L**3 + (Q(699)/2 + 180*s5)*L**2*lam + (-93 - 48*s5)*L*lam**2,
 4: (Q(1)/6 + s5/2)*Li(4, 1/(2*phi**2)) + (-Q(163)/48 + Q(13)/16*s5)*Li(4, phi**-2) + (Q(29)/30 - Q(11)/10*s5)*Li(4, Q(1)/2)
    + (Q(11)/2 + s5/2)*Li(4, 2/phi**2) + (Q(11)/3 + s5)*Li(4, phi/2) + (-Q(1)/8 + s5/72)*Li(4, phi**-6) + (1 + Q(7)/9*s5)*Li(4, -phi**-3)
    + (Q(52)/5 + Q(18)/5*s5)*L*z3 + (-Q(27)/5 - Q(8)/5*s5)*lam*z3 + (-9 - 4*s5)*L*Li(3, 1/phi) + (4 + 2*s5)*lam*Li(3, 1/phi)
    + (-Q(191)/2700 - Q(47)/8100*s5)*mp.pi**4 + (-Q(155)/72 - Q(13)/24*s5)*pi2*L**2 + (Q(113)/45 + s5/3)*pi2*L*lam + (-Q(2)/3 - s5/30)*pi2*lam**2
    + (Q(1939)/144 + Q(59)/48*s5)*L**4 + (-Q(385)/18 + s5/6)*L**3*lam + (Q(133)/12 - Q(3)/4*s5)*L**2*lam**2 + (-Q(43)/18 + s5/6)*L*lam**3 + (Q(17)/60 + s5/60)*lam**4})
for tag, (c, forms) in cf.items():
    W = split(c); M3 = sum(W.values())
    ref = mp.mpf(open(C9 + f"{tag}_M3_hpl_1150.txt").read().split()[0])
    print(f"{tag}: c0 = {mp.nstr(c, 12)}  M3(planner) = {mp.nstr(M3, 30)}  |M3 - COR9 file| = {mp.nstr(abs(M3 - ref), 3)}  weights {sorted(W)}")
    for w in (2, 3, 4):
        print(f"   W{w} = {mp.nstr(W[w], 22)}   |planner split - COR9 closed form| = {mp.nstr(abs(W[w] - forms[w]), 3)}")
