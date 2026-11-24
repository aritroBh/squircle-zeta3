
import mpmath as mp, json, time
mp.mp.dps = 60
t0 = time.time()
out = {}

# ---------- Part A: the slope law ----------
def phi(s, th):
    c = mp.cos(th)
    return 1 - (s + c)/mp.sqrt(1 + 2*c*s + s**2)
def S_of_U(u, th):
    C = mp.cos(th/2)**2; Z = mp.sin(th/2)**2
    return C/u - mp.cos(th) - Z*u
def Lambda(nu, th):
    z = mp.sin(th/2)**2
    return 2*z*mp.hyper([-nu, 1+nu], [2], z)

# A1: mass = 1-cos theta
d1 = max(abs(mp.quad(lambda s: phi(s,th), [0,mp.inf]) - (1-mp.cos(th)))
         for th in [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, mp.mpf('2.7')])
# A2: Mellin identity E[S^w] (via explicit law) = (pi w/sin pi w) Lambda_w/(1-c)
ths = [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, mp.mpf('2.7')]
ws  = [mp.mpf('-0.62'), mp.mpf('0.3'), mp.mpf('0.88')]
d2 = mp.mpf(0)
for th in ths:
    c = mp.cos(th)
    for w in ws:
        lhs = mp.quad(lambda u: S_of_U(u, th)**w, [0,1])
        rhs = (mp.pi*w/mp.sin(mp.pi*w))*Lambda(w, th)/(1-c)
        d2 = max(d2, abs(lhs-rhs))
out['A_mass'] = mp.nstr(d1,3); out['A_mellin'] = mp.nstr(d2,3)

# A3: tail identity P(S>s) (explicit U root) = defective tail / (1-c)
d3 = mp.mpf(0)
for th in ths:
    c = mp.cos(th); Z = mp.sin(th/2)**2; C = mp.cos(th/2)**2
    for s0 in [mp.mpf('0.4'), 1, 5]:
        u_star = (-(s0+c) + mp.sqrt((s0+c)**2 + 4*Z*C))/(2*Z)
        tail = mp.quad(lambda s: phi(s,th), [s0, mp.inf])/(1-c)
        d3 = max(d3, abs(tail - u_star))
out['A_tail'] = mp.nstr(d3,3)

# A4: geometric mean E log S_theta = a1 = -1 - cot^2 log cos^2
d4 = mp.mpf(0)
for th in ths:
    a1 = -1 - (1/mp.tan(th/2)**2)*mp.log(mp.cos(th/2)**2)
    num = mp.quad(lambda u: mp.log(S_of_U(u,th)), [0,1])
    d4 = max(d4, abs(a1-num))
out['A_geomean'] = mp.nstr(d4,3)
print('Part A (slope law):', out, f'{time.time()-t0:.1f}s')
json.dump(out, open('partA.json','w'))
