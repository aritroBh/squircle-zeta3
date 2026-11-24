code = r'''
import mpmath as mp, pickle, time
import sympy as sp
mp.mp.dps = 620
t0=time.time()
D1 = '../DOD1_dodecagon/'
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
S = pickle.load(open(D1+'S_assembly.pkl','rb'))['S']
cf = pickle.load(open('atom_coefs_pq.pkl','rb'))
r3 = sp.sqrt(3)
# build full substitution dict
X = atoms['X']
sub = {}
for i in X:
    sub[sp.symbols('L2x%d'%i)] = atoms['L2'][i]
    sub[sp.symbols('L3x%d'%i)] = atoms['L3'][i]
    sub[sp.symbols('Ax%d'%i)] = atoms['Ax'][i]
for k in atoms['T2']:
    sub[sp.symbols('T2_%d_%d'%k)] = atoms['T2'][k]
    sub[sp.symbols('T3p_%d_%d'%k)] = atoms['T3p'][k]
for k in atoms['Lam']:
    sub[sp.symbols('Lam%d_%d'%k)] = atoms['Lam'][k]
for k in atoms['T3a']:
    sub[sp.symbols('T3_%d_%d%d'%k)] = atoms['T3a'][k]
sub[sp.symbols('l2')] = atoms['l2']; sub[sp.symbols('l3')] = atoms['l3']
sub[sp.symbols('lu1')] = atoms['lu']; sub[sp.symbols('z2')] = atoms['z2']; sub[sp.symbols('z3')] = atoms['z3']
print("subs dict size:", len(sub), "t=%ds"%(time.time()-t0))
S_full = sp.N(S.subs(sub, simultaneous=True), 620)
print("S_full evaluated t=%ds"%(time.time()-t0))
# zero out Lam and T3a
sub0 = dict(sub)
for k in atoms['Lam']: sub0[sp.symbols('Lam%d_%d'%k)] = mp.mpf(0)
for k in atoms['T3a']: sub0[sp.symbols('T3_%d_%d%d'%k)] = mp.mpf(0)
S_nod2 = sp.N(S.subs(sub0, simultaneous=True), 620)
print("S_nod2 evaluated t=%ds"%(time.time()-t0))
Xi_from_S = S_full - S_nod2
print("Xi(S_full - S_nod2) =", mp.nstr(Xi_from_S, 20))
# compare with Xi_num from extracted coefficients
Lam = atoms['Lam']; T3a = atoms['T3a']
XiD = mp.mpf(0); XiD_p = mp.mpf(0); XiD_q = mp.mpf(0)
for (fam,ij) in cf:
    if fam=='Lam': v = Lam[ij]
    elif fam=='T3a': v = T3a[ij]
    else: continue
    p,q = cf[(fam,ij)]
    XiD += (mp.mpf(str(p))+mp.mpf(str(q))*mp.sqrt(3))*v
    XiD_p += mp.mpf(str(p))*v; XiD_q += mp.mpf(str(q))*v
print("XiD (extracted)       =", mp.nstr(XiD, 20))
print("|Xi_from_S - XiD| =", mp.nstr(abs(Xi_from_S - XiD),4))
print("XiD_p =", mp.nstr(XiD_p,20), " XiD_q =", mp.nstr(XiD_q,20))
print("|XiD - (XiD_p+sqrt3 XiD_q)| =", mp.nstr(abs(XiD-(XiD_p+mp.sqrt(3)*XiD_q)),4))
# also check S_pq_620 stored consistency
spq = pickle.load(open(D1+'S_pq_620.pkl','rb'))
print("stored S_p =", mp.nstr(spq['S_p'],20), " S_q =", mp.nstr(spq['S_q'],20))
print("S_full - (S_p + sqrt3 S_q) =", mp.nstr(abs(S_full - (spq['S_p']+mp.sqrt(3)*spq['S_q'])),4))
pickle.dump({'Xi_from_S':Xi_from_S,'XiD':XiD,'XiD_p':XiD_p,'XiD_q':XiD_q,'S_full':S_full,'S_nod2':S_nod2},
            open('Xi_verified_620.pkl','wb'))
print("elapsed", time.time()-t0)
'''
open('step_11.py','w').write(code)
exec(code)
