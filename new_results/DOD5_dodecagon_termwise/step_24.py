import pickle, mpmath as mp
mp.mp.dps = 620
base = "<project>/research/08_theory/"
D1, D2, D3 = base+"DOD1_dodecagon/", base+"DOD2_dodecagon_finish/", base+"DOD3_dodecagon_symbol/"
ac = pickle.load(open(D2+"atom_coefs_pq.pkl","rb"))
xv = pickle.load(open(D2+"Xi_verified_620.pkl","rb"))
Rv = pickle.load(open(D3+"R_all12_vals620.pkl","rb"))['Rq']
atoms = pickle.load(open(D1+"atoms620.pkl","rb"))
g = pickle.load(open(D2+"gdswaps620.pkl","rb"))
r3 = mp.sqrt(3)

# what Xi36 should be: Xi - sum c^T3a * R
cR = mp.mpf(0)
for (i,j,k) in sorted(Rv):
    p,q = ac[('T3a',(i,j,k))]
    cR += (mp.mpf(p.numerator)/p.denominator + r3*mp.mpf(q.numerator)/q.denominator)*Rv[(i,j,k)]
Xi36_target = mp.mpf(xv['Xi_from_S']) - cR
print("Xi36 target (= Xi - sum cR):", mp.nstr(Xi36_target, 40))

# what DOD4 step_60 computed: X0v + r3*X1v from atoms/gdswaps with coefficients
IDX=[1,2,5,6]
Xv = mp.mpf(0)
for (kind, ij),(p,q) in ac.items():
    c = mp.mpf(p.numerator)/p.denominator + r3*mp.mpf(q.numerator)/q.denominator
    if kind=='Lam':
        i,j = ij
        Xv += c*atoms['Lam'][(i,j)]
    else:
        i,j,k = ij
        xi = atoms['X'][i]
        Xv += c/xi*(g['SwG'][(i,j,k)] - g['SwD'][(i,j,k)])
print("DOD4-assembled Xi36 (atoms/gdswaps):", mp.nstr(Xv, 40))
print("diff:", mp.nstr(abs(Xv - Xi36_target), 3))
