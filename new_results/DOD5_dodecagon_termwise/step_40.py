import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
from dod4_lib import *
d1p = pickle.load(open('dod5_task1.pkl','rb'))
TERMS, vals = d1p['TERMS'], d1p['vals']
mp.mp.dps = 620
atoms = pickle.load(open(D1+'atoms620.pkl','rb'))
ac = pickle.load(open(D2+"atom_coefs_pq.pkl","rb"))
xv = pickle.load(open(D2+"Xi_verified_620.pkl","rb"))
rv = pickle.load(open(D3+"R_all12_vals620.pkl","rb"))
r3 = mp.sqrt(3)

# Xi36 with accurate swaps
Xi36_acc = mp.mpf(0)
for (kind, ij),(p,q) in ac.items():
    c = mp.mpf(p.numerator)/p.denominator + r3*mp.mpf(q.numerator)/q.denominator
    if kind=='Lam':
        Xi36_acc += c*vals[('Lam',ij)]
    else:
        i,j,k = ij
        Xi36_acc += c/atoms['X'][i]*(vals[('SwG',ij)] - vals[('SwD',ij)])
# Xi - sum c*Rv(sym)
cR = mp.mpf(0)
for k in rv['Rv']:
    p,q = ac[('T3a',k)]
    cR += (mp.mpf(p.numerator)/p.denominator + r3*mp.mpf(q.numerator)/q.denominator)*rv['Rv'][k]
target = mp.mpf(xv['Xi_from_S']) - cR
print("Xi36_acc =", mp.nstr(Xi36_acc, 40))
print("target   =", mp.nstr(target, 40))
print("|Xi36_acc - (Xi - sum cR_sym)| =", mp.nstr(abs(Xi36_acc - target), 3))

# Also: Rq_accurate vs Rv sym, already done. Save accurate values + R
R_acc = {}
TRIPLES = [(1,2,5),(1,2,6),(1,5,6),(2,1,5),(2,1,6),(2,5,6),(5,1,2),(5,1,6),(5,2,6),(6,1,2),(6,1,5),(6,2,5)]
for (i,j,k) in TRIPLES:
    R_acc[(i,j,k)] = atoms['T3a'][(i,j,k)] - (vals[('SwG',(i,j,k))]-vals[('SwD',(i,j,k))])/atoms['X'][i]
pickle.dump({'Xi36_acc':Xi36_acc,'R_acc':R_acc}, open('dod5_xi36_acc.pkl','wb'))
print("saved dod5_xi36_acc.pkl")
# and confirm Xi36_acc + sum cR = Xi to 1e-618
tot = Xi36_acc + sum((mp.mpf(ac[('T3a',k)][0].numerator)/ac[('T3a',k)][0].denominator + r3*mp.mpf(ac[('T3a',k)][1].numerator)/ac[('T3a',k)][1].denominator)*R_acc[k] for k in TRIPLES)
print("|Xi36_acc + sum cR_acc - Xi| =", mp.nstr(abs(tot - mp.mpf(xv['Xi_from_S'])),3))
