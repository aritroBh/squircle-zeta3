import os, mpmath as mp, pickle
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 830
v = pickle.load(open('atoms830.pkl','rb'))
lf,l2,z2,z3,phi = v['lf'],v['l2'],v['z2'],v['z3'],v['phi']
L3x4,L3x1,L3xf = v['L3_x4'],v['L3_x1'],v['L3_xf']
names = ['z3','z2*lf','z2*l2','lf^3','lf^2*l2','lf*l2^2','l2^3','L3(-1/phi)','L3(1/phi)']
B = [z3, z2*lf, z2*l2, lf**3, lf**2*l2, lf*l2**2, l2**3, L3x1, L3xf]
rels={}
for tname in ['A_x1','A_x4','Lam','swap']:
    rel = mp.pslq([v[tname]]+B, tol=mp.mpf('1e-790'), maxcoeff=10**14, maxsteps=10000)
    if rel is None:
        print(tname, "NO relation"); rels[tname]=None; continue
    resid = rel[0]*v[tname] + sum(rel[i+1]*B[i] for i in range(len(B)))
    rels[tname]=rel
    print(tname, rel)
    print("   resid:", mp.nstr(abs(resid),8))
pickle.dump(rels, open('rels_lll2.pkl','wb'))
with open('progress.log','a') as f:
    f.write(f"step: pslq vs independent 9-const basis: A_x1 {rels['A_x1']}, A_x4 {rels['A_x4']}, Lam {rels['Lam']}, swap {rels['swap']}\n")