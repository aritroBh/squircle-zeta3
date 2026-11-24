import sympy as sp, pickle, os, sys, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
import dec_engine as E
r5 = sp.sqrt(5)
mons = pickle.load(open('mons_dec2.pkl','rb'))
Jloc = pickle.load(open('Jloc_dec2.pkl','rb'))
Jtab = pickle.load(open('Jtab_dec2.pkl','rb'))
MIX = pickle.load(open('MIX20_proved.pkl','rb'))
def JLOC(n,a,b):
    return sp.expand(Jloc[(n,a,b)].subs(sp.pi**2, 6*E.z2).subs(sp.zeta(3), E.z3))
def K1v(n,a,b,c,side):
    return sp.expand(Jtab[(n,a,b,c)].subs(E.sub1 if side==1 else E.sub4, simultaneous=True))
def INT(n,a,b,c1,c4):
    if c1==0 and c4==0: return JLOC(n,a,b)
    if c4==0: return K1v(n,a,b,c1,1)
    if c1==0: return K1v(n,a,b,c4,4)
    return MIX[(n,a,b,c1,c4)]
t0=time.time()
l2, lf = E.l2, E.lf
S = sp.Integer(0)
for (n_, abc, letters), coef in mons.items():
    a,b,c1,c4 = abc
    S += coef * l2**letters[0] * lf**letters[1] * INT(n_,a,b,c1,c4)
S = sp.expand(S)
print("assembly done %.0fs; terms:"%(time.time()-t0), len(sp.Add.make_args(S)))
# --- substitutions (all PROVED) ---
subs_red = {
  E.L2_1: -2*E.z2/5 + lf**2/2,          # golden Li2 (Euler+Landen+duplication, certified)
  E.L2_4:  2*E.z2/5 - lf**2,
  E.L2f:   3*E.z2/5 - lf**2,
  E.AA_1:  4*E.z3/5 - 2*lf*E.z2/5 + 2*E.L3_1,        # DEC1 orbit proof
  E.AA_4:  E.z3 - 3*lf*E.z2/5 + 4*E.L3_1 + 3*E.L3f,  # DEC1 orbit proof
  E.Lam:   2*E.z3/5 + 4*lf*E.z2/5 - lf**3 + 3*E.L3_1 + E.L3f,   # PRV1 Thm D
  E.swap:  3*E.z3/5 - 7*lf*E.z2/5 + 2*lf**3 + E.L3_1 + E.L3f,   # PRV1 Thm D
  E.L3_4:  4*E.L3_1 + 4*E.L3f,           # duplication Ra (exact series)
}
S3 = sp.expand(S.subs(subs_red, simultaneous=True))
pickle.dump({'S':S, 'S3':S3}, open('S_assembly_dec2.pkl','wb'))
print("reduced; basis vars:", set(S3.free_symbols))
