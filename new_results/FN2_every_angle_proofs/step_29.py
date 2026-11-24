import mpmath as mp, time
mp.mp.dps = 60
t0 = time.time()

def phi_st(s, c):
    # stable form: phi_theta(s) = (1-c^2)/( sqrt(1+2cs+s^2) * (sqrt(1+2cs+s^2)+s+c) )
    r = mp.sqrt(1 + 2*c*s + s**2)
    return (1-c**2)/(r*(r+s+c))

def S12(z):
    return mp.zeta(3) - mp.polylog(3, 1-z) + mp.log(1-z)*mp.polylog(2, 1-z) \
           + mp.mpf('0.5')*mp.log(1-z)**2*mp.log(z)
def cumulants_closed(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    return (1+u, 1+mp.pi**2/3-2*mp.polylog(2,z)-u**2, -12*w*S12(z)+2*(1+u**3))

def moments_logsub(th, rmax=4):
    c = mp.cos(th); out = []
    for r in range(1, rmax+1):
        i1 = mp.quad(lambda u: phi_st(mp.exp(-u),c)*mp.exp(-u)*u**r, [0, mp.inf])
        i2 = mp.quad(lambda u: phi_st(mp.exp(u),c)*mp.exp(u)*(-u)**r, [0, mp.inf])
        out.append((i1+i2)/(1-c))
    return out
def raw_to_cumulants(m):
    m1,m2,m3,m4 = m
    return (m1, m2-m1**2, m3-3*m2*m1+2*m1**3, m4-4*m3*m1-3*m2**2+12*m2*m1**2-6*m1**4)

worst = 0; k4tab = {}
print("== Lemma C (stable form): cumulants of log(1/S_theta) ==")
for th in [mp.pi/6, mp.pi/4, mp.pi/3, mp.pi/2, 2*mp.pi/3, 3*mp.pi/4, 5*mp.pi/6]:
    kn = raw_to_cumulants(moments_logsub(th))
    kc = cumulants_closed(th)
    res = [abs(kn[j]-kc[j]) for j in range(3)]
    worst = max(worst, *res)
    k4tab[mp.nstr(th,6)] = mp.nstr(kn[3],12)
    print(f"th={mp.nstr(th,8)}: |k1|={mp.nstr(res[0],2)} |k2|={mp.nstr(res[1],2)} |k3|={mp.nstr(res[2],2)}  k4={mp.nstr(kn[3],10)}")

def mellin_num(w, th):
    c = mp.cos(th)
    i1 = mp.quad(lambda u: phi_st(mp.exp(-u),c)*mp.exp(-u)*mp.exp(-w*u), [0, mp.inf])
    i2 = mp.quad(lambda u: phi_st(mp.exp(u),c)*mp.exp(u)*mp.exp(w*u), [0, mp.inf])
    return (i1+i2)/(1-c)
def mellin_legendre(w, th):
    c = mp.cos(th); z = mp.sin(th/2)**2
    return mp.pi*w/mp.sin(mp.pi*w)*(2*z*mp.hyper([-w,1+w],[2],z))/(1-c)
print("\nMellin check (stable):")
wm = 0
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3, 5*mp.pi/6]:
    for w in [mp.mpf('-0.95'), mp.mpf('-0.5'), mp.mpf('0.5'), mp.mpf('0.95'), mp.mpf('0.3')]:
        r = abs(mellin_num(w,th)-mellin_legendre(w,th)); wm = max(wm,r)
        print(f"  th={mp.nstr(th,6)} w={mp.nstr(w,5)}: |res|={mp.nstr(r,2)}")
print(f"\nworst cumulant resid = {mp.nstr(worst,2)}; worst Mellin resid = {mp.nstr(wm,2)}")
print("elapsed", time.time()-t0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} Lemma C stable: worst cumulant resid {mp.nstr(worst,2)} over 7 angles, Mellin worst {mp.nstr(wm,2)}  [NUMERICAL]\n")
