
import mpmath as mp
mp.mp.dps = 60
def S_of_U(u, th):
    C = mp.cos(th/2)**2; Z = mp.sin(th/2)**2
    return C/u - mp.cos(th) - Z*u
def Lambda(nu, th):
    z = mp.sin(th/2)**2
    return 2*z*mp.hyper([-nu, 1+nu], [2], z)
ths = [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, mp.mpf('2.7')]
d3 = mp.mpf(0)
for th in ths:
    c = mp.cos(th)
    for w in [mp.mpf('-0.62'), mp.mpf('0.3'), mp.mpf('0.88')]:
        # u = e^{-t}: integrand S(e^{-t})^w e^{-t}, t in [0, inf), decays like e^{-(1-w)t}
        lhs = mp.quad(lambda t: S_of_U(mp.exp(-t), th)**w * mp.exp(-t), [0, 3, 10, mp.inf], maxdegree=9)
        rhs = (mp.pi*w/mp.sin(mp.pi*w))*Lambda(w, th)/(1-c)
        d3 = max(d3, abs(lhs-rhs))
print('A3 Mellin max diff (u=e^-t):', mp.nstr(d3,3))
import json
out = json.load(open('partA.json')); out['A3_mellin'] = mp.nstr(d3,3); json.dump(out, open('partA.json','w'))
print(out)
