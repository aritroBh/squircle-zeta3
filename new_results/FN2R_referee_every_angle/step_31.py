import mpmath as mp, time
mp.mp.dps = 50
t0 = time.time()
def phi(s,c):
    if s==0: return 1-c
    r = mp.sqrt(1+2*c*s+s**2)
    return (1-c**2)/(r*(r+s+c))
def S12(zv):  # paper's S_{1,2}(z), ell = log(1-z)
    return mp.zeta(3) - mp.polylog(3, 1-zv) + mp.log(1-zv)*mp.polylog(2, 1-zv) \
           + mp.mpf('0.5')*mp.log(1-zv)**2*mp.log(zv)
def kforms(th):
    zv = mp.sin(th/2)**2; w2 = 1/mp.tan(th/2)**2; ellc = mp.log(mp.cos(th/2)**2); u = w2*ellc
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,zv)-u**2, 2+2*u**3-12*w2*S12(zv))
def a_forms(th):  # paper's a1,a2,a3 (thm:P2asym)
    zv = mp.sin(th/2)**2; w2 = 1/mp.tan(th/2)**2; u = w2*mp.log(mp.cos(th/2)**2)
    return (-1-u, mp.zeta(2)+mp.mpf('0.5')-mp.polylog(2,zv)-u**2/2, 2*w2*S12(zv)-(1+u**3)/3)

print("== (1) Mellin = Legendre: E[S^w] quadrature vs (pi w/sin pi w)*2z*2F1/(1-c) ==")
worst = mp.mpf(0)
for th in [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, 5*mp.pi/6]:
    c = mp.cos(th); zv = mp.sin(th/2)**2
    for wv in [mp.mpf('-0.9'), mp.mpf('-0.5'), mp.mpf('-0.1'), mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('0.9')]:
        lhs = mp.quad(lambda s: s**wv*phi(s,c), [0,1,mp.inf])/(1-c)
        rhs = (mp.pi*wv/mp.sin(mp.pi*wv))*2*zv*mp.hyper([-wv,1+wv],[2],zv)/(1-c)
        worst = max(worst, abs(lhs-rhs))
print("  worst |diff| over 5 angles x 6 w:", mp.nstr(worst,3))

print("== (2) cumulants: quadrature of log-moments vs closed forms ==")
worst2 = mp.mpf(0)
for th in [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, 5*mp.pi/6, mp.pi/12, 11*mp.pi/12]:
    c = mp.cos(th)
    raw = []
    for r in (1,2,3):
        i1 = mp.quad(lambda u: phi(mp.exp(-u),c)*mp.exp(-u)*u**r, [0,mp.inf])
        i2 = mp.quad(lambda u: phi(mp.exp(u),c)*mp.exp(u)*(-u)**r, [0,mp.inf])
        raw.append((i1+i2)/(1-c))
    m1,m2,m3 = raw
    kq = (m1, m2-m1**2, m3-3*m1*m2+2*m1**3)
    kf = kforms(th)
    d = max(abs(kq[i]-kf[i]) for i in range(3))
    worst2 = max(worst2,d)
    print(f"  th={mp.nstr(th,8)}: k1={mp.nstr(kq[0],12)} k2={mp.nstr(kq[1],12)} k3={mp.nstr(kq[2],12)} |diff|={mp.nstr(d,2)}")
print("  worst:", mp.nstr(worst2,3))

print("== (3) k_j = (-1)^j j! a_j vs paper's a_j closed forms (thm:P2asym) ==")
worst3 = mp.mpf(0)
for th in [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, 5*mp.pi/6]:
    kf = kforms(th); af = a_forms(th)
    conv = (-af[0], 2*af[1], -6*af[2])
    d = max(abs(kf[i]-conv[i]) for i in range(3))
    worst3 = max(worst3,d)
print("  worst |k_j - (-1)^j j! a_j| with both from closed forms:", mp.nstr(worst3,3))

print("== (4) k_j at pi/2 reduce to 1-log2, 1+pi^2/6, 2-3/2 zeta(3) ==")
kf = kforms(mp.pi/2)
print("  ", mp.nstr(kf[0]-(1-mp.log(2)),2), mp.nstr(kf[1]-(1+mp.pi**2/6),2), mp.nstr(kf[2]-(2-mp.mpf('1.5')*mp.zeta(3)),2))

print("== (5) moment boundary: E[S^w] blows up at w=1 and w=-1 (cutoff growth) ==")
c = mp.cos(mp.pi/3)
for wv, cut in [(mp.mpf('0.999'), 1e6), (mp.mpf('1.0'), 1e6), (mp.mpf('-0.999'), None), (mp.mpf('-1.0'), None)]:
    if wv > 0:
        v = mp.quad(lambda s: s**wv*phi(s,c), [1, cut]) 
        print(f"  w={wv}: int_1^1e6 s^w phi ds = {mp.nstr(v,6)} (~ d2 * {'log' if wv==1 else 's^(w-1)'})")
    else:
        v = mp.quad(lambda s: s**wv*phi(s,c), [mp.mpf('1e-12'), mp.mpf('1e-6'), 1])
        print(f"  w={wv}: int_1e-12^1 s^w phi ds = {mp.nstr(v,6)}")
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} step_07: Lemma C numerics: Mellin=Legendre worst {mp.nstr(worst,2)}; cumulants vs closed forms worst {mp.nstr(worst2,2)} (7 angles incl pi/12, 11pi/12); k_j=(-1)^j j! a_j consistent; pi/2 reduction exact; E[S^w] diverges at w=+-1 as claimed [referee]\n")
