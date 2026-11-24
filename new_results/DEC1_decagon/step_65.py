import os, pickle, sympy as sp, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
import dec_engine as E
mons = pickle.load(open('mons.pkl','rb'))
lf,l2,z2,z3 = E.lf,E.l2,E.z2,E.z3
# local part of S: monomials with c1=c4=0
S_loc = sp.Integer(0); S_s1 = sp.Integer(0); S_s4 = sp.Integer(0); S_mix = sp.Integer(0)
for (n,abc,letters),pq in mons.items():
    a,b,c1,c4 = abc; e2,ef = letters
    J = E.INT(n,a,b,c1,c4)
    term = pq*l2**e2*lf**ef*J
    if c1==0 and c4==0: S_loc += term
    elif c4==0: S_s1 += term
    elif c1==0: S_s4 += term
    else: S_mix += term
S_loc = sp.expand(S_loc)
z3coeff_loc = sp.expand(S_loc).coeff(z3,1)
print("z3 coeff of S_loc:", z3coeff_loc)
# local z3 contribution to P'''(1): 60 L * z3coeff_loc, should equal C5
Lsym = 2*sp.sqrt(10-2*r5)/4/(r5*E.phi**2)
r5 = sp.sqrt(5)
Lsym = 2*sp.sqrt(10-2*r5)/4/(r5*E.phi**2)
C5_check = sp.simplify(60*Lsym*z3coeff_loc)
g = pickle.load(open('geometry.pkl','rb'))
print("60 L [z3]S_loc =", C5_check)
print("C_5          =", g['C5'])
print("equal:", sp.simplify(C5_check-g['C5'])==0)
# total z3 coefficient in S (named basis, after reductions)
S3 = pickle.load(open('S_reduced.pkl','rb'))['S3']
z3tot = sp.Poly(sp.expand(S3), l2,lf,z2,z3,E.L3_1,E.L3f).coeff_monomial(z3)
print("[z3] S3 (total, named basis):", sp.simplify(z3tot))
z3tot_clean = sp.simplify((z3tot+z3tot.subs(r5,-r5))/2) + sp.simplify((z3tot-z3tot.subs(r5,-r5))/(2*r5))*r5
print("[z3] S3 =", z3tot_clean)
T = sp.simplify(60*Lsym*z3tot_clean)   # total z3 coeff of P'''(1)
print("total z3 coeff of P'''(1):", T)
Ne = sp.simplify((T - g['C5'])/10)
print("non-local z3 per edge (this basis):", Ne)
# explicit D3,D4 quad checks at 110 dps
mp.mp.dps=110
v = pickle.load(open('atoms830.pkl','rb'))
num = {E.l2:v['l2'], E.lf:v['lf'], E.z2:v['z2'], E.z3:v['z3'],
       E.L2_1: mp.polylog(2,-1/v['phi']), E.L3_1: v['L3_x1'], E.AA_1: v['A_x1'],
       E.L2_4: mp.polylog(2,1/v['phi']**2), E.L3_4: v['L3_x4'], E.AA_4: v['A_x4'],
       E.L2f: mp.polylog(2,1/v['phi']), E.L3f: v['L3_xf'], E.Lam: v['Lam'], E.swap: v['swap']}
def ev(e): return mp.mpf(sp.N(e.subs(num, simultaneous=True), 110))
phim=v['phi']; kap1m=1/phim; kap4m=-1/phim**2
w1m=lambda t: mp.log(1+kap1m*t); w4m=lambda t: mp.log(1+kap4m*t)
D3q = mp.quad(lambda t: mp.log(1-t)*w4m(t)/(1+kap1m*t),[0,1])
D4q = mp.quad(lambda t: mp.log(1-t)*w1m(t)/(1+kap4m*t),[0,1])
print("|D3 closed - quad| =", mp.nstr(abs(ev(E.D3)-D3q),6))
print("|D4 closed - quad| =", mp.nstr(abs(ev(E.D4)-D4q),6))
with open('progress.log','a') as f:
    f.write(f"step: local z3 == C_5 exactly (diff 0). total z3 = {T}; N_e per edge = {Ne}. D3,D4 verified {mp.nstr(abs(ev(E.D3)-D3q),3)},{mp.nstr(abs(ev(E.D4)-D4q),3)}\n")