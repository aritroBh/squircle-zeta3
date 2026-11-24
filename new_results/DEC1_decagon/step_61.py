import os, pickle, sympy as sp, mpmath as mp, time
os.chdir('research/08_theory/DEC1_decagon')
import dec_engine as E
mons = pickle.load(open('mons.pkl','rb'))
r5 = sp.sqrt(5)
t0=time.time()
S = sp.Integer(0)
for (n,abc,letters),pq in mons.items():
    a,b,c1,c4 = abc
    e2,ef = letters
    J = E.INT(n,a,b,c1,c4)
    S += pq * E.l2**e2 * E.lf**ef * J
S = sp.expand(S)
print("assembly time %.1fs, terms:"%(time.time()-t0), len(sp.Add.make_args(S)))
pickle.dump({'S':S}, open('S_assembly.pkl','wb'))
# numeric evaluation
mp.mp.dps = 110
v = pickle.load(open('atoms830.pkl','rb'))
phim=v['phi']; lfm=v['lf']; l2m=v['l2']; z2m=v['z2']; z3m=v['z3']
num = {E.l2:l2m, E.lf:lfm, E.z2:z2m, E.z3:z3m,
       E.L2_1: mp.polylog(2,-1/phim), E.L3_1: v['L3_x1'], E.AA_1: v['A_x1'],
       E.L2_4: mp.polylog(2,1/phim**2), E.L3_4: v['L3_x4'], E.AA_4: v['A_x4'],
       E.L2f: mp.polylog(2,1/phim), E.L3f: v['L3_xf'], E.Lam: v['Lam'], E.swap: v['swap']}
Sn = mp.mpf(sp.N(S.subs(num, simultaneous=True), 110))
s36 = mp.sqrt(10-2*mp.sqrt(5))/4
Lm = 2*s36/(mp.sqrt(5)*phim**2)
P3_asm = 60*Lm*Sn
gold = pickle.load(open('P3_660.pkl','rb'))['P3']
print("P3 assembled:", mp.nstr(P3_asm,50))
print("gold        :", mp.nstr(gold,50))
print("diff:", mp.nstr(abs(P3_asm-gold),6))
with open('progress.log','a') as f:
    f.write(f"step: full assembly S ({len(sp.Add.make_args(S))} terms) -> P3_asm; |P3_asm - 640d gold| = {mp.nstr(abs(P3_asm-gold),4)} (110d eval)\n")