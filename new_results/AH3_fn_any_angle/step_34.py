
import mpmath as mp, time
mp.mp.dps = 60
def S_of_U(u, th):
    C = mp.cos(th/2)**2; Z = mp.sin(th/2)**2
    return C/u - mp.cos(th) - Z*u
def Lambda(nu, th):
    z = mp.sin(th/2)**2
    return 2*z*mp.hyper([-nu, 1+nu], [2], z)
ths = [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, mp.mpf('2.7')]
d3 = mp.mpf(0); worst=None
for th in ths:
    c = mp.cos(th)
    for w in [mp.mpf('-0.62'), mp.mpf('0.3'), mp.mpf('0.88')]:
        lhs = mp.quad(lambda u: S_of_U(u, th)**w, [0, mp.mpf('1e-5'), mp.mpf('1e-3'), mp.mpf('0.1'), mp.mpf('0.5'), 1], maxdegree=10)
        rhs = (mp.pi*w/mp.sin(mp.pi*w))*Lambda(w, th)/(1-c)
        if abs(lhs-rhs) > d3: d3 = abs(lhs-rhs); worst=(mp.nstr(th,6), w, mp.nstr(lhs,20), mp.nstr(rhs,20))
print('A3 Mellin max diff:', mp.nstr(d3,3), '\nworst:', worst)
