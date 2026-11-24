import os, pickle, sympy as sp, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
import dec_engine as E
sp_t = sp.symbols('t')
mp.mp.dps = 110
v = pickle.load(open('atoms830.pkl','rb'))
phim=v['phi']; lfm=v['lf']; l2m=v['l2']; z2m=v['z2']; z3m=v['z3']
num = {E.l2:l2m, E.lf:lfm, E.z2:z2m, E.z3:z3m,
       E.L2_1: mp.polylog(2,-1/phim), E.L3_1: v['L3_x1'], E.AA_1: v['A_x1'],
       E.L2_4: mp.polylog(2,1/phim**2), E.L3_4: v['L3_x4'], E.AA_4: v['A_x4'],
       E.L2f: mp.polylog(2,1/phim), E.L3f: v['L3_xf'], E.Lam: v['Lam'], E.swap: v['swap']}
def ev(e): return mp.mpf(sp.N(e.subs(num, simultaneous=True), 110))
kap1m=1/phim; kap4m=-1/phim**2
w1m=lambda t: mp.log(1+kap1m*t); w4m=lambda t: mp.log(1+kap4m*t)
worst=mp.mpf(0); nb=0
# Jloc
for (n,a,b) in sorted(E.Jloc.keys()):
    f=lambda t,n=n,a=a,b=b: t**n*(mp.log(t))**a*(mp.log(1-t))**b
    d = abs(ev(E.JLOC(n,a,b))-mp.quad(f,[0,1])); worst=max(worst,d); nb+=1
print("Jloc verified:", nb, "worst", mp.nstr(worst,4))
# Jtab at x1 (side1) and x4 (side4)
worst1=mp.mpf(0); worst4=mp.mpf(0)
for (n,a,b,c) in sorted(E.Jtab.keys()):
    f1=lambda t,n=n,a=a,b=b,c=c: t**n*(mp.log(t))**a*(mp.log(1-t))**b*w1m(t)**c
    f4=lambda t,n=n,a=a,b=b,c=c: t**n*(mp.log(t))**a*(mp.log(1-t))**b*w4m(t)**c
    worst1=max(worst1,abs(ev(E.K1(n,a,b,c,1))-mp.quad(f1,[0,1])))
    worst4=max(worst4,abs(ev(E.K1(n,a,b,c,4))-mp.quad(f4,[0,1])))
print("Jtab x1 worst:", mp.nstr(worst1,4), " Jtab x4 worst:", mp.nstr(worst4,4))
with open('progress.log','a') as f:
    f.write(f"step: Jloc (40) worst {mp.nstr(worst,4)}; Jtab at x1,x4 (40 each) worst {mp.nstr(worst1,4)}, {mp.nstr(worst4,4)} (110 dps quad)\n")