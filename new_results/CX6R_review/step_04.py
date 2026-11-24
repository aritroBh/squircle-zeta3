import mpmath as mp
mp.mp.dps = 50
w = mp.expjpi(mp.mpf(1)/4)
a = mp.pi/(2*mp.sqrt(2))
def H(s): return mp.sqrt(mp.pi)/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf(3)/2)
def H_int(s):  # independent: definition integral
    return mp.quad(lambda u: mp.sqrt(1-mp.exp(-2*u))*mp.exp(-s*u), [0,1,5,mp.inf])
def K(s): return 1/(2*mp.cos(mp.pi*s/2))
def E(s): return mp.expj(mp.pi*s/2)
def C(z): return mp.exp(z)*(1+mp.exp(-z)-mp.sqrt(1+mp.exp(-2*z)))
print("CX6R step_04: sec.2 tests, 50 dps")
# (a) H beta identity
dH = max(abs(H(s)-H_int(s)) for s in [mp.mpf('0.5'), 1, 3, mp.mpf('2.3')+1j, mp.mpf('0.5')+2*w*3])
print("max |H_gamma - H_integral| at 5 pts:", mp.nstr(dH,6))
# (b) bounds (4): |H| <= 1/x and <= sqrt(pi/2) x^{-3/2}
worst1 = worst2 = mp.mpf(0)
for i in range(1,200):
    x = mp.mpf('0.05')*i
    s = x + mp.mpf('3.7')*1j*x  # generic complex s with Re s = x
    hv = abs(H(s))
    worst1 = max(worst1, hv*x)
    worst2 = max(worst2, hv*x**mp.mpf('1.5'))
print("max |H|*Re(s) (must be <=1):", mp.nstr(worst1,10))
print("max |H|*Re(s)^1.5 (must be <= sqrt(pi/2)=%s):"%mp.nstr(mp.sqrt(mp.pi/2),10), mp.nstr(worst2,10))
# (c) identity (5): C(z) = 1 - (2/pi) Im int_Gamma H K e^{-sz} ds
def gamma_int(z):
    f = lambda r: H(mp.mpf('0.5')+r*w)*K(mp.mpf('0.5')+r*w)*mp.exp(-(mp.mpf('0.5')+r*w)*z)*w
    return mp.quad(f, [0, 2, 8, 30, mp.inf])
for z in [mp.mpf('0'), mp.mpf('0.3'), mp.mpf('1'), mp.mpf('2.5'), mp.mpf('10')]:
    lhs = C(z)
    rhs = 1 - 2/mp.pi*mp.im(gamma_int(z))
    # also series (6): sum_n (-1)^n H(2n+1) e^{-(2n+1)z}
    ser = mp.nsum(lambda n: (-1)**n*H(2*n+1)*mp.exp(-(2*n+1)*z), [0, mp.inf])
    print(f"z={mp.nstr(z,6)}  C(z)={mp.nstr(lhs,25)}")
    print(f"   |C-(1-2/pi Im int)|={mp.nstr(abs(lhs-rhs),6)}   |(2/pi)series-(1-C)|={mp.nstr(abs(2/mp.pi*ser-(1-lhs)),6)}")
# (d) residue check: Res_{s=2n+1} K = (-1)^{n+1}/pi
for n in range(4):
    s0 = 2*n+1
    res = mp.residue(K, s0)
    print(f"Res K at {s0}: {mp.nstr(res,20)}  vs (-1)^{n+1}/pi = {mp.nstr((-1)**(n+1)/mp.pi,20)}")
# (e) bounds (13) on Gamma, fine grid
m1 = m2 = m3 = mp.mpf(0)
for i in range(2000):
    r = mp.mpf('0.01')*i
    s = mp.mpf('0.5')+r*w
    m1 = max(m1, abs(H(s)))
    m2 = max(m2, abs(K(s))*mp.exp(a*r))
    m3 = max(m3, abs(K(s)-E(s))*mp.exp(3*a*r))
print("max |H| on Gamma grid (bound 2):", mp.nstr(m1,10))
print("max |K| e^{ar} (bound 2):", mp.nstr(m2,10))
print("max |K-E| e^{3ar} (bound 2):", mp.nstr(m3,10))
# true minimum of |1+e^{i pi s}| on Gamma (proof needs >= 1/2)
mn = min(abs(1+mp.exp(mp.j*mp.pi*(mp.mpf('0.5')+mp.mpf('0.001')*i*w))) for i in range(20000))
print("min |1+e^{i pi s}| on Gamma (proof uses >= 1/2):", mp.nstr(mn,15))
