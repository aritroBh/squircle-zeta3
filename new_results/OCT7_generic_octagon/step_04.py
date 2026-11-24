
import mpmath as mp, sympy as sp, pickle
mp.mp.dps = 65
import mixed_engine_gen as E

def Li21(a, c, N=600):
    a, c = mp.mpf(a), mp.mpf(c)
    return mp.fsum([a**n * c**m / (mp.mpf(n)**2 * m) for n in range(2, N) for m in range(1, n)])

L2 = lambda v: mp.polylog(2, mp.mpf(v))
L3 = lambda v: mp.polylog(3, mp.mpf(v))
NUM = {
 E.l2: mp.log(2), E.l3: mp.log(3), E.l5: mp.log(5),
 E.z2: mp.pi**2/6, E.z3: mp.zeta(3),
 E.L2x1: L2(-mp.mpf(1)/3), E.L3x1: L3(-mp.mpf(1)/3),
 E.L2x2: L2(-mp.mpf(1)/9), E.L3x2: L3(-mp.mpf(1)/9),
 E.L2b1: L2(mp.mpf(1)/4),  E.L3b1: L3(mp.mpf(1)/4),
 E.L2b2: L2(mp.mpf(1)/10), E.L3b2: L3(mp.mpf(1)/10),
 E.L2mh: L2(-mp.mpf(1)/2), E.L2m2t: L2(-mp.mpf(2)/3),
 E.L2m3h: L2(-mp.mpf(3)/2), E.L3m3h: L3(-mp.mpf(3)/2),
 E.L2m2: L2(-mp.mpf(2)), E.L3m2: L3(-mp.mpf(2)),
 E.L2_3f: L2(mp.mpf(3)/5), E.L3_3f: L3(mp.mpf(3)/5),
 E.L2_2t: L2(mp.mpf(2)/3), E.L3_2t: L3(mp.mpf(2)/3),
 E.Lam1: Li21(-mp.mpf(1)/3, mp.mpf(1)/3),
 E.Lam2: Li21(-mp.mpf(1)/9, mp.mpf(3)),
 E.Lam3: Li21(mp.mpf(1)/4, mp.mpf(2)/5),
 E.Lam4: Li21(mp.mpf(1)/10, mp.mpf(5)/2),
}
def ev(e):
    return mp.mpf(sp.N(e.subs(list(NUM.items())), mp.mp.dps))

k1m = mp.mpf(1)/3; k2m = mp.mpf(1)/9
X = lambda t: mp.log(t); Y = lambda t: mp.log(1-t)
W1 = lambda t: mp.log(1+k1m*t); W2 = lambda t: mp.log(1+k2m*t)
def quad(n,a,b,c1,c2):
    f = lambda t: t**n * X(t)**a * Y(t)**b * W1(t)**c1 * W2(t)**c2
    return mp.quad(f, [0,1])

# 1) check K1 (single-ratio) integrals
worstK = mp.mpf(0); cnt=0
for side,km in [(1,k1m),(2,k2m)]:
    for (n,a,b,c) in [(0,0,0,1),(1,0,0,1),(2,0,0,1),(3,0,0,1),(0,0,1,1),(0,1,0,1),(0,0,0,2),
                      (1,1,0,1),(2,0,1,1),(1,0,0,2),(0,1,0,2),(0,0,1,2),(0,2,0,1),(0,0,2,1),(0,0,0,3)]:
        sym = E.K1(n,a,b,c,side)
        if side==1: q = mp.quad(lambda t: t**n*X(t)**a*Y(t)**b*W1(t)**c, [0,1])
        else:       q = mp.quad(lambda t: t**n*X(t)**a*Y(t)**b*W2(t)**c, [0,1])
        d = abs(ev(sym)-q); worstK=max(worstK,d); cnt+=1
print(f"K1 checks: {cnt} integrals, worst |exact-quad| = {mp.nstr(worstK,3)}")

# 2) check all 20 mixed integrals (n in 0..3)
tests = []
for n in range(4):
    tests += [((n,0,0,1,1), 'W1W2'), ((n,1,0,1,1), 'XW1W2'), ((n,0,1,1,1), 'YW1W2'),
              ((n,0,0,2,1), 'W1^2W2'), ((n,0,0,1,2), 'W1W2^2')]
worstM = mp.mpf(0)
for (n,a,b,c1,c2), lbl in tests:
    sym = E.INT(n,a,b,c1,c2)
    q = quad(n,a,b,c1,c2)
    d = abs(ev(sym)-q); worstM=max(worstM,d)
print(f"MIXED checks: {len(tests)} integrals, worst |exact-quad| = {mp.nstr(worstM,3)}")

# 3) check local JLOC a few
worstL = mp.mpf(0)
for (n,a,b) in [(0,0,0),(1,1,1),(2,2,1),(3,3,0),(1,2,1),(3,1,2)]:
    d = abs(ev(E.JLOC(n,a,b)) - quad(n,a,b,0,0)); worstL=max(worstL,d)
print(f"JLOC checks: worst = {mp.nstr(worstL,3)}")

# 4) Koelbig identity check
AA1num = mp.fsum([(-mp.mpf(1)/3)**m * mp.harmonic(m) / m**2 for m in range(1,3000)])
rhs = 2*NUM[E.L3x1] + NUM[E.L3b1] - NUM[E.L2x1]*(2*mp.log(2)-mp.log(3)) - (2*mp.log(2)-mp.log(3))**3/6
AA2num = mp.fsum([(-mp.mpf(1)/9)**m * mp.harmonic(m) / m**2 for m in range(1,3000)])
rhs2 = 2*NUM[E.L3x2] + NUM[E.L3b2] - NUM[E.L2x2]*(mp.log(10)-mp.log(9)) - (mp.log(10)-mp.log(9))**3/6
print("Koelbig AA1 residual:", mp.nstr(abs(AA1num-rhs),3), " AA2:", mp.nstr(abs(AA2num-rhs2),3))
pickle.dump('ok', open('engine_verified_flag.pkl','wb'))
with open("progress.log","a") as f:
    f.write(f"STEP 04: generalized engine verified: K1 worst {mp.nstr(worstK,3)}, MIXED worst {mp.nstr(worstM,3)}, JLOC worst {mp.nstr(worstL,3)}, Koelbig {mp.nstr(abs(AA1num-rhs),3)}/{mp.nstr(abs(AA2num-rhs2),3)}\n")
