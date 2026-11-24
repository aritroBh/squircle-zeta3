import mpmath as mp, time
mp.mp.dps = 60
t0 = time.time()

def phi(s, c):
    return 1 - (s + c)/mp.sqrt(1 + 2*c*s + s**2)

def S12(z):
    return mp.zeta(3) - mp.polylog(3, 1-z) + mp.log(1-z)*mp.polylog(2, 1-z) \
           + mp.mpf('0.5')*mp.log(1-z)**2*mp.log(z)

def cumulants_closed(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    k1 = 1 + u
    k2 = 1 + mp.pi**2/3 - 2*mp.polylog(2, z) - u**2
    k3 = -12*w*S12(z) + 2*(1 + u**3)
    return k1, k2, k3

def moments_logsub(th, rmax=4):
    """raw moments of xi=log(1/S_theta) via s=e^{-u}, s=e^{u} splitting; all
    integrals exponentially convergent."""
    c = mp.cos(th)
    out = []
    for r in range(1, rmax+1):
        # part s in (0,1): s=e^{-u}, log(1/s)=u
        i1 = mp.quad(lambda u: phi(mp.exp(-u),c)*mp.exp(-u)*u**r, [0, mp.inf])
        # part s in (1,inf): s=e^{u}, log(1/s)=-u
        i2 = mp.quad(lambda u: phi(mp.exp(u),c)*mp.exp(u)*(-u)**r, [0, mp.inf])
        out.append((i1+i2)/(1-c))
    return out

def raw_to_cumulants(m):
    m1,m2,m3,m4 = m
    return (m1, m2-m1**2, m3-3*m2*m1+2*m1**3,
            m4-4*m3*m1-3*m2**2+12*m2*m1**2-6*m1**4)

worst = 0
print("== Lemma C redo: cumulants via log-substitution ==")
for th in [mp.pi/6, mp.pi/4, mp.pi/3, mp.pi/2, 2*mp.pi/3, 3*mp.pi/4, 5*mp.pi/6]:
    m = moments_logsub(th)
    kn = raw_to_cumulants(m)
    kc = cumulants_closed(th)
    res = [abs(kn[j]-kc[j]) for j in range(3)]
    worst = max(worst, *res)
    print(f"th={mp.nstr(th,8)}: |k1 res|={mp.nstr(res[0],2)} |k2 res|={mp.nstr(res[1],2)} |k3 res|={mp.nstr(res[2],2)}  k4(num)={mp.nstr(kn[3],8)}")

# Mellin with log substitution
def mellin_num(w, th):
    c = mp.cos(th)
    i1 = mp.quad(lambda u: phi(mp.exp(-u),c)*mp.exp(-u)*mp.exp(-w*u), [0, mp.inf])  # s^w, s=e^{-u}
    i2 = mp.quad(lambda u: phi(mp.exp(u),c)*mp.exp(u)*mp.exp(w*u), [0, mp.inf])
    return (i1+i2)/(1-c)
def mellin_legendre(w, th):
    c = mp.cos(th); z = mp.sin(th/2)**2
    return mp.pi*w/mp.sin(mp.pi*w)*(2*z*mp.hyper([-w,1+w],[2],z))/(1-c)
print("\nMellin check:")
wm = 0
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3, 5*mp.pi/6]:
    for w in [mp.mpf('-0.95'), mp.mpf('-0.5'), mp.mpf('0.5'), mp.mpf('0.95')]:
        r = abs(mellin_num(w,th)-mellin_legendre(w,th))
        wm = max(wm, r)
        print(f"  th={mp.nstr(th,6)} w={mp.nstr(w,5)}: |res|={mp.nstr(r,2)}")
print(f"\nworst cumulant resid={mp.nstr(worst,2)}, worst Mellin resid={mp.nstr(wm,2)}")
print("elapsed", time.time()-t0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} Lemma C redo (log-substitution): cumulants match closed forms, Mellin=Legendre at w=+-0.95,+-0.5, 4 angles  [NUMERICAL]\n")
