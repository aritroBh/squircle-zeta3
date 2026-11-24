import mpmath as mp, sympy as sp
mp.mp.dps = 65
import mixed_engine_gen as E

def Li21(a, c, N=600):
    a, c = mp.mpf(a), mp.mpf(c)
    return mp.fsum([a**n * c**m / (mp.mpf(n)**2 * m) for n in range(2, N) for m in range(1, n)])
L2 = lambda v: mp.polylog(2, mp.mpf(v)); L3 = lambda v: mp.polylog(3, mp.mpf(v))
NUM = {
 E.l2: mp.log(2), E.l3: mp.log(3), E.l5: mp.log(5), E.z2: mp.pi**2/6, E.z3: mp.zeta(3),
 E.L2x1: L2(-mp.mpf(1)/3), E.L3x1: L3(-mp.mpf(1)/3), E.L2x2: L2(-mp.mpf(1)/9), E.L3x2: L3(-mp.mpf(1)/9),
 E.L2b1: L2(mp.mpf(1)/4), E.L3b1: L3(mp.mpf(1)/4), E.L2b2: L2(mp.mpf(1)/10), E.L3b2: L3(mp.mpf(1)/10),
 E.L2mh: L2(-mp.mpf(1)/2), E.L2m2t: L2(-mp.mpf(2)/3), E.L2m3h: L2(-mp.mpf(3)/2), E.L3m3h: L3(-mp.mpf(3)/2),
 E.L2m2: L2(-mp.mpf(2)), E.L3m2: L3(-mp.mpf(2)), E.L2_3f: L2(mp.mpf(3)/5), E.L3_3f: L3(mp.mpf(3)/5),
 E.L2_2t: L2(mp.mpf(2)/3), E.L3_2t: L3(mp.mpf(2)/3),
 E.Lam1: Li21(-mp.mpf(1)/3, mp.mpf(1)/3), E.Lam2: Li21(-mp.mpf(1)/9, mp.mpf(3)),
 E.Lam3: Li21(mp.mpf(1)/4, mp.mpf(2)/5), E.Lam4: Li21(mp.mpf(1)/10, mp.mpf(5)/2),
}
def ev(e): return mp.mpf(sp.N(e.subs(list(NUM.items())), mp.mp.dps))
k1m = mp.mpf(1)/3; k2m = mp.mpf(1)/9
X = lambda t: mp.log(t); Y = lambda t: mp.log(1-t)
W1 = lambda t: mp.log(1+k1m*t); W2 = lambda t: mp.log(1+k2m*t)
def quad(n,a,b,c1,c2):
    return mp.quad(lambda t: t**n * X(t)**a * Y(t)**b * W1(t)**c1 * W2(t)**c2, [0,1])

for n in range(2):
    for (a,b,c1,c2,lbl) in [(0,0,1,1,'MixedW'),(1,0,1,1,'M1'),(0,1,1,1,'M2'),(0,0,2,1,'M3'),(0,0,1,2,'M4')]:
        sym = E.INT(n,a,b,c1,c2); q = quad(n,a,b,c1,c2)
        print(f"n={n} {lbl}: exact={mp.nstr(ev(sym),12)} quad={mp.nstr(q,12)} diff={mp.nstr(abs(ev(sym)-q),3)}")
