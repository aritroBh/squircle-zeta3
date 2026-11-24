import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
import mpmath as mp
mp.mp.dps = 620
import sympy as sp
t0=time.time()
out = pickle.load(open('dod5_atom_coefs_all.pkl','rb'))
sd = pickle.load(open(base+"DOD1_dodecagon/S_atom_decomp.pkl","rb"))
atoms = pickle.load(open(base+"DOD1_dodecagon/atoms620.pkl","rb"))
xv = pickle.load(open(base+"DOD2_dodecagon_finish/Xi_verified_620.pkl","rb"))
r3 = mp.sqrt(3)
l2,l3,lu1,z2,z3 = sp.symbols('l2 l3 lu1 z2 z3')

def fpq(pq):
    p,q = pq
    return mp.mpf(str(p)) + r3*mp.mpf(str(q))
# S_singles_part value
ss = sd['S_singles_part']
subv = {l2: mp.log(2), l3: mp.log(3), lu1: mp.log(2-r3), z2: mp.pi**2/6, z3: mp.zeta(3)}
ss_val = mp.mpf(sp.N(ss.subs(subv), 620))
print("S_singles_part value:", mp.nstr(ss_val, 30))

S_d1 = ss_val
for k in out:
    fam, kk = k
    if fam in ('Lam','T3a'): continue
    if fam=='T2': v = atoms['T2'][kk]
    elif fam=='T3p': v = atoms['T3p'][kk]
    elif fam=='Ax': v = atoms['Ax'][kk]
    elif fam=='L2x': v = atoms['L2'][kk]
    elif fam=='L3x': v = atoms['L3'][kk]
    cd = out[k]
    cval = fpq(cd.get('1',(0,0)))
    for ls in ('l2','l3','lu1'):
        if ls in cd: cval += fpq(cd[ls])*subv[sp.Symbol(ls)]
    S_d1 += cval*v
print("S_d1 assembled =", mp.nstr(S_d1, 40))
print("S_nod2         =", mp.nstr(mp.mpf(xv['S_nod2']), 40))
print("|diff| =", mp.nstr(abs(S_d1 - mp.mpf(xv['S_nod2'])),3), " t=%.0fs"%(time.time()-t0))
pickle.dump({'ss_val':ss_val,'S_d1':S_d1}, open('dod5_sd1_check.pkl','wb'))
