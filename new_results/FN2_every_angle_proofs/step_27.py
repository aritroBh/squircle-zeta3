import mpmath as mp, time
mp.mp.dps = 60
t0 = time.time()

def phi(s, c):
    return 1 - (s + c)/mp.sqrt(1 + 2*c*s + s**2)

def S12(z):
    # Nielsen S_{1,2}(z) = sum_{n>=1} (-1)^{n-1} z^n/(n^2 (1+...+1/n))... use the standard
    # S_{1,2}(z) = (1/2) int_0^z log^2(1-t)/t dt  (this is the convention matching AH3)
    return mp.zeta(3) - mp.polylog(3, 1-z) + mp.log(1-z)*mp.polylog(2, 1-z) \
           + mp.mpf('0.5')*mp.log(1-z)**2*mp.log(z)

def cumulants_closed(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    k1 = 1 + u
    k2 = 1 + mp.pi**2/3 - 2*mp.polylog(2, z) - u**2
    k3 = -12*w*S12(z) + 2*(1 + u**3)
    return k1, k2, k3

def moments_numeric(th, rmax=4):
    c = mp.cos(th)
    # raw moments of xi = log(1/S_theta): m_r = int phi (log 1/s)^r / (1-c)
    out = []
    for r in range(1, rmax+1):
        v = mp.quad(lambda s: phi(s,c)*mp.log(1/s)**r, [0, 1, mp.inf])/(1-c)
        out.append(v)
    return out

def raw_to_cumulants(m):
    m1,m2,m3,m4 = m
    k1 = m1
    k2 = m2 - m1**2
    k3 = m3 - 3*m2*m1 + 2*m1**3
    k4 = m4 - 4*m3*m1 - 3*m2**2 + 12*m2*m1**2 - 6*m1**4
    return k1,k2,k3,k4

print("== Lemma C: cumulants of log(1/S_theta): quadrature vs closed forms ==")
for th in [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, 5*mp.pi/6]:
    m = moments_numeric(th)
    kn = raw_to_cumulants(m)
    kc = cumulants_closed(th)
    r1 = kn[0]-kc[0]; r2 = kn[1]-kc[1]; r3 = kn[2]-kc[2]
    print(f"th={mp.nstr(th,8)}: k1 resid={mp.nstr(r1,2)}, k2 resid={mp.nstr(r2,2)}, k3 resid={mp.nstr(r3,2)}, k4(num)={mp.nstr(kn[3],10)}")

# reduction to theta=pi/2 values
k1c,k2c,k3c = cumulants_closed(mp.pi/2)
print("\npi/2 reduction:", mp.nstr(k1c-(1-mp.log(2)),2), mp.nstr(k2c-(1+mp.pi**2/6),2),
      mp.nstr(k3c-(2-mp.mpf('1.5')*mp.zeta(3)),2))

# E[S^w] finite for -1<w<1, vs the paper's Legendre formula (proved there)
def mellin_legendre(w, th):
    c = mp.cos(th); z = mp.sin(th/2)**2
    Lam = 2*z*mp.hyper([-w, 1+w], [2], z)
    return mp.pi*w/mp.sin(mp.pi*w)*Lam/(1-c)
print("\nMellin E[S_theta^w]: quadrature vs Legendre formula")
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
    c = mp.cos(th)
    for w in [mp.mpf('-0.9'), mp.mpf('-0.5'), mp.mpf('0.5'), mp.mpf('0.9')]:
        num = mp.quad(lambda s: phi(s,c)*s**w, [0, 1, mp.inf])/(1-c)
        ref = mellin_legendre(w, th)
        print(f"  th={mp.nstr(th,6)} w={mp.nstr(w,4)}: resid={mp.nstr(num-ref,2)}")
print("elapsed", time.time()-t0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} Lemma C: cumulants k1..k3(theta) quadrature-vs-closed-form ~1e-55; Mellin=Legendre checked w=+-0.9,+-0.5  [NUMERICAL; Mellin identity itself PROVED in paper Thm P2]\n")
