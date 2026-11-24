import sympy as sp, pickle, os, sys
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
r5 = sp.sqrt(5)
tt = sp.Symbol('t')
exact = pickle.load(open('E3n_sym_pq.pkl','rb'))
L0,L1,L2,L3,L4 = sp.symbols('L0 L1 L2 L3 L4')
lt,l1t,w1s,w4s,l2,lf = sp.symbols('lt l1t w1 w4 l2 lf')
Ls=[L0,L1,L2,L3,L4]
E = sp.Integer(0)
for key,row in exact.items():
    mon = sp.Integer(1)
    for i,e in enumerate(key): mon *= Ls[i]**e
    for nn in range(4):
        p,q = row[nn]
        if p==0 and q==0: continue
        E += (p+q*r5)*tt**nn*mon
E = sp.expand(E)
sub = {L0: lt-l2-2*lf, L1: l1t-l2-2*lf, L2: w1s-l2-2*lf, L3: w4s-l2-lf, L4: -l2-lf}
E2 = sp.expand(E.xreplace(sub))
mons = {}
for term in sp.Add.make_args(E2):
    d = term.as_powers_dict()
    n_ = int(d.get(tt,0))
    abc = (int(d.get(lt,0)), int(d.get(l1t,0)), int(d.get(w1s,0)), int(d.get(w4s,0)))
    letters = (int(d.get(l2,0)), int(d.get(lf,0)))
    rest = sp.simplify(term/(tt**n_ * lt**abc[0]*l1t**abc[1]*w1s**abc[2]*w4s**abc[3] * l2**letters[0]*lf**letters[1]))
    key = (n_, abc, letters)
    mons[key] = sp.expand(mons.get(key,0)+rest)
print("my mons:", len(mons))
dec1 = pickle.load(open('<project>/research/08_theory/DEC1_decagon/mons.pkl','rb'))
print("DEC1 mons:", len(dec1))
mism = 0
for k in dec1:
    if k not in mons:
        print("missing", k); mism += 1; continue
    d = sp.simplify(mons[k] - dec1[k])
    if d != 0:
        mism += 1
        if mism <= 5: print("MISMATCH", k, d, "mine:", mons[k], "dec1:", dec1[k])
extra = set(mons)-set(dec1)
print("keys only in mine:", extra)
print("coefficient mismatches:", mism)
pickle.dump(mons, open('mons_dec2.pkl','wb'))
