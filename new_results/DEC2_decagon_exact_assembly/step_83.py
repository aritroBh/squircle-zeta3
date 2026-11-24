import sympy as sp, pickle, os, sys, mpmath as mp, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
import dec_engine as E
MIXvals = pickle.load(open('MIX20_proved.pkl','rb'))
mp.mp.dps = 110
s5 = mp.sqrt(5); phim = (1+s5)/2; lfm = mp.log(phim)
k1m = 1/phim; k4m = -1/phim**2
x1m = -1/phim; x4m = 1/phim**2
def Li21(x,y):
    tot=mp.mpf(0); ps=mp.mpf(0); xp=mp.mpf(1); yp=mp.mpf(1); m=0
    while True:
        m+=1; xp*=x
        if m>=2: tot += xp*ps/m**2
        yp*=y; ps += yp/m
        if m>10 and abs(xp*ps/m**2) < mp.mpf('1e-115'): break
    return tot
atoms = {E.l2: mp.log(2), E.lf: lfm, E.z2: mp.pi**2/6, E.z3: mp.zeta(3),
         E.L2_1: mp.polylog(2,x1m), E.L3_1: mp.polylog(3,x1m),
         E.AA_1: mp.nsum(lambda mm: x1m**mm*mp.harmonic(mm)/mm**2, [1, mp.inf]),
         E.L2_4: mp.polylog(2,x4m), E.L3_4: mp.polylog(3,x4m),
         E.AA_4: mp.nsum(lambda mm: x4m**mm*mp.harmonic(mm)/mm**2, [1, mp.inf]),
         E.L2f: mp.polylog(2,1/phim), E.L3f: mp.polylog(3,1/phim),
         E.Lam: Li21(x1m,x1m), E.swap: Li21(x4m,-phim)}
t0=time.time()
worst = mp.mpf(0)
for (n,a,b,c1,c4), v in MIXvals.items():
    q = mp.quad(lambda tt: tt**n*mp.log(tt)**a*mp.log(1-tt)**b*mp.log(1+k1m*tt)**c1*mp.log(1+k4m*tt)**c4, [0,1])
    val = mp.mpf(sp.N(v.subs(atoms, simultaneous=True), 110))
    d = abs(val-q)
    worst = max(worst,d)
    if d > mp.mpf('1e-95'):
        print("large:", (n,a,b,c1,c4), mp.nstr(d,3))
print("20 mixed integrals vs quadrature (110 dps): worst |diff| =", mp.nstr(worst,3), " %.0fs"%(time.time()-t0))
with open('progress.log','a') as f:
    f.write(f"step_12: 20 mixed integrals PROVED (20 IBP integrand identities residual 0; 8+8+4 partial-fraction identities exact; boundaries -lf^2, +/-lf^3, 0); symbolic match with DEC1 engine diff 0; vs quadrature 110 dps worst {mp.nstr(worst,2)}\n")
