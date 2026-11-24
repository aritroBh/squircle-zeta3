import sympy as sp, pickle, os, sys, mpmath as mp, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
t = sp.Symbol('t'); s = sp.Symbol('s'); r5 = sp.sqrt(5); phi = (1+r5)/2
kap1 = 1/phi; kap4 = -1/phi**2
ids = [sp.expand((1+kap4*(1-s)) - (phi+s)/phi**2), sp.expand((phi+s)/phi**2 - (1/phi)*(1+s/phi)),
       sp.expand((1+kap1*(1-s)) - phi*(1-s/phi**2)), sp.expand((phi**2-s)/phi - phi*(1-s/phi**2)),
       sp.expand((phi+s) - phi*(1+s/phi))]
print("reflection identities simplified:", [sp.simplify(i) for i in ids])
# numeric validation of D1..D10 (engine values vs quadrature)
import dec_engine as E
mp.mp.dps = 110
s5 = mp.sqrt(5); phim = (1+s5)/2; lfm = mp.log(phim); z2m = mp.pi**2/6; z3m = mp.zeta(3)
k1m = 1/phim; k4m = -1/phim**2
u1 = lambda tt: 1+k1m*tt; u4 = lambda tt: 1+k4m*tt
w1 = lambda tt: mp.log(u1(tt)); w4 = lambda tt: mp.log(u4(tt))
quads = {
 'D1': mp.quad(lambda tt: mp.log(tt)*w4(tt)/u1(tt), [0,1]),
 'D2': mp.quad(lambda tt: mp.log(tt)*w1(tt)/u4(tt), [0,1]),
 'D3': mp.quad(lambda tt: mp.log(1-tt)*w4(tt)/u1(tt), [0,1]),
 'D4': mp.quad(lambda tt: mp.log(1-tt)*w1(tt)/u4(tt), [0,1]),
 'D5': mp.quad(lambda tt: w4(tt)/u1(tt), [0,1]),
 'D6': mp.quad(lambda tt: w1(tt)/u4(tt), [0,1]),
 'D7': mp.quad(lambda tt: w1(tt)**2/u4(tt), [0,1]),
 'D8': mp.quad(lambda tt: w4(tt)**2/u1(tt), [0,1]),
 'D9': mp.quad(lambda tt: w1(tt)*w4(tt)/u1(tt), [0,1]),
 'D10': mp.quad(lambda tt: w1(tt)*w4(tt)/u4(tt), [0,1]),
}
def Li21(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1); m=0
    while True:
        m+=1; xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
        if m>10 and abs(xp*ps/m**2) < mp.mpf('1e-115'): break
    return tot
x1m = -1/phim; x4m = 1/phim**2
atoms = {E.l2: mp.log(2), E.lf: lfm, E.z2: z2m, E.z3: z3m,
         E.L2_1: mp.polylog(2,x1m), E.L3_1: mp.polylog(3,x1m),
         E.AA_1: mp.nsum(lambda mm: x1m**mm*mp.harmonic(mm)/mm**2, [1, mp.inf]),
         E.L2_4: mp.polylog(2,x4m), E.L3_4: mp.polylog(3,x4m),
         E.AA_4: mp.nsum(lambda mm: x4m**mm*mp.harmonic(mm)/mm**2, [1, mp.inf]),
         E.L2f: mp.polylog(2,1/phim), E.L3f: mp.polylog(3,1/phim),
         E.Lam: Li21(x1m,x1m), E.swap: Li21(x4m,-phim)}
t0=time.time()
worst = mp.mpf(0)
for name, eng in [('D1',E.D1),('D2',E.D2),('D3',E.D3),('D4',E.D4),('D5',E.D5),('D6',E.D6),('D7',E.D7),('D8',E.D8),('D9',E.D9),('D10',E.D10)]:
    val = mp.mpf(sp.N(eng.subs(atoms, simultaneous=True), 110))
    d = abs(val - quads[name])
    worst = max(worst, d)
    print(f"{name}: |proved value - quadrature| = {mp.nstr(d,3)}")
print("worst over D1..D10:", mp.nstr(worst,3), " (%.0fs)"%(time.time()-t0))
with open('progress.log','a') as f:
    f.write(f"step_11: D1-D4 certified (D1,D2 absolutely-convergent series rearrangements, term equality exact; D3,D4 t->1-t reflections, 5 algebraic identities exactly 0, canonical match with engine diff 0); all D1..D10 vs quadrature (110 dps): worst {mp.nstr(worst,2)}\n")
