
import mpmath as mp, json
mp.mp.dps = 90
def S_of_U(u, th):
    C = mp.cos(th/2)**2; Z = mp.sin(th/2)**2
    return C/u - mp.cos(th) - Z*u
def Lambda_int(nu, th):
    return mp.quad(lambda t: mp.legenp(nu, 0, t), [mp.cos(th), 1])
ths = [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, mp.mpf('2.7')]
d3 = mp.mpf(0)
for th in ths:
    c = mp.cos(th)
    for w in [mp.mpf('-0.62'), mp.mpf('0.3'), mp.mpf('0.88')]:
        lhs = mp.quad(lambda t: S_of_U(mp.exp(-t), th)**w * mp.exp(-t), [0, 3, 10, 25, mp.inf], maxdegree=10)
        rhs = (mp.pi*w/mp.sin(mp.pi*w))*Lambda_int(w, th)/(1-c)
        d3 = max(d3, abs(lhs-rhs))
print('A3 Mellin max diff at dps=90:', mp.nstr(d3,3))
out = json.load(open('partA.json')); out['A3_mellin_dps90'] = mp.nstr(d3,3); json.dump(out, open('partA.json','w'))
