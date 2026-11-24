
import mpmath as mp, json, time
mp.mp.dps = 60
t0=time.time()

def phi(s, th):
    c = mp.cos(th)
    return 1 - (s + c)/mp.sqrt(1 + 2*c*s + s**2)
def S_of_U(u, th):
    C = mp.cos(th/2)**2; Z = mp.sin(th/2)**2
    return C/u - mp.cos(th) - Z*u
def Lambda(nu, th):
    z = mp.sin(th/2)**2
    return 2*z*mp.hyper([-nu, 1+nu], [2], z)
def rho(s, c):
    return mp.sqrt(1 + 2*c*s + s**2) - s - c

ths = [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, mp.mpf('2.7')]
out = {}

# A1 mass: tail rho(s) at s=0 should be 1-c ; rho(infty)=0
d1 = max(abs(rho(mp.mpf(0), mp.cos(th)) - (1-mp.cos(th))) for th in ths)
# A2 tail identity: rho_c(s) = (1-c) u*(s)  (pure algebra, 60 digits)
d2 = mp.mpf(0)
for th in ths:
    c = mp.cos(th); Z = mp.sin(th/2)**2; C = mp.cos(th/2)**2
    for s0 in [mp.mpf('0.01'), mp.mpf('0.4'), 1, 5, 400]:
        u_star = (-(s0+c) + mp.sqrt((s0+c)**2 + 4*Z*C))/(2*Z)
        d2 = max(d2, abs(rho(s0,c) - (1-c)*u_star))
# A3 Mellin identity (split quadrature, higher degree)
d3 = mp.mpf(0)
for th in ths:
    c = mp.cos(th)
    for w in [mp.mpf('-0.62'), mp.mpf('0.3'), mp.mpf('0.88')]:
        lhs = mp.quad(lambda u: S_of_U(u, th)**w, [0, mp.mpf('0.02'), mp.mpf('0.5'), 1], maxdegree=10)
        rhs = (mp.pi*w/mp.sin(mp.pi*w))*Lambda(w, th)/(1-c)
        d3 = max(d3, abs(lhs-rhs))
# A4 geometric mean
d4 = mp.mpf(0)
for th in ths:
    a1 = -1 - (1/mp.tan(th/2)**2)*mp.log(mp.cos(th/2)**2)
    num = mp.quad(lambda u: mp.log(S_of_U(u,th)), [0, mp.mpf('1e-8'), mp.mpf('0.5'), 1], maxdegree=10)
    d4 = max(d4, abs(a1-num))
out = {'A1_mass': mp.nstr(d1,3), 'A2_tail': mp.nstr(d2,3), 'A3_mellin': mp.nstr(d3,3), 'A4_geomean': mp.nstr(d4,3)}
print('Part A:', out, f'{time.time()-t0:.1f}s')
json.dump(out, open('partA.json','w'))

# ---------- Part B: cumulants, 60 digits ----------
def S12(z):
    ell = mp.log(1-z)
    return mp.zeta(3) - mp.polylog(3, 1-z) + ell*mp.polylog(2, 1-z) + mp.mpf('0.5')*ell**2*mp.log(z)
def cumulants_closed(th):
    z = mp.sin(th/2)**2; w = 1/mp.tan(th/2)**2; ell = mp.log(mp.cos(th/2)**2); u = w*ell
    k1 = 1 + u
    k2 = 1 + mp.pi**2/3 - 2*mp.polylog(2, z) - u**2
    k3 = -12*w*S12(z) + 2*(1 + u**3)
    return k1, k2, k3
def cumulants_numeric(th, J=3):
    f = lambda t: mp.log(mp.quad(lambda u: S_of_U(u, th)**(-t), [0, mp.mpf('1e-6'), mp.mpf('0.3'), 1], maxdegree=10))
    return [mp.diff(f, mp.mpf(0), j) for j in range(1, J+1)]
d5 = mp.mpf(0)
for th in ths:
    kc = cumulants_closed(th); kn = cumulants_numeric(th)
    for j in range(3): d5 = max(d5, abs(kc[j]-kn[j]))
print('Part B (cumulants k1..k3 closed vs numeric):', mp.nstr(d5,3))
out['B_cumulants'] = mp.nstr(d5,3)
json.dump(out, open('partA.json','w'))
