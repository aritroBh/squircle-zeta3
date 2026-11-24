import sys, pickle, json, random
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
import evsetup
from evsetup import cc, h_k, GAL, a2mp
import mpmath as mp
mp.mp.dps=130
evsetup.setup(700,1000,800)
BASE='<project>'
SQ2=mp.sqrt(2)

# --- load test relations ---
rogers=json.load(open(BASE+'/research/08_theory/A3f_level8_relations/weight2_new_relation.json'))['row']
rogers=[(tuple(w),c) for w,c in rogers]  # integer coeffs
five=pickle.load(open(BASE+'/research/08_theory/A3f4_glanois_basis/five_relations.pickle','rb'))
fw3=five['freewords']; R=five['relations']
Ts=pickle.load(open(BASE+'/research/08_theory/A3h_weight4_target/w4_systems_T.pkl','rb'))
T3=Ts['T3']
print('T3 rows:',len(T3),'type of row:',type(T3[0]), T3[0][0])
std_small=[(lab,row) for lab,row in T3 if lab[0]=='std' and len(row)<=6][:3]
for lab,row in std_small: print('std row',lab,'support',len(row),'sample',row[:3])

def evalrel(terms,k):
    """terms: list of (word, (a,b)) meaning (b + a*sqrt2)*h(w). Returns value at embedding k."""
    s=GAL[k]
    tot=mp.mpc(0)
    for w,(a,b) in terms:
        c = mp.mpf(str(b)) + s*SQ2*mp.mpf(str(a))
        tot += c*a2mp(h_k(w,k),120)
    return tot

print()
print('=== Rogers five-term (weight 2) at k=1,3,5,7 ===')
terms=[(w,(0,c)) for w,c in rogers]
for k in (1,3,5,7):
    v=evalrel(terms,k); print(f'k={k}: |v|={mp.nstr(abs(v),4)}')
print('=== R1..R5 (weight 3) at 4 embeddings ===')
for i,r in enumerate(R):
    terms=[(tuple(fw3[j]),(0,r[j])) for j in range(32)]
    vs=[abs(evalrel(terms,k)) for k in (1,3,5,7)]
    print(f'R{i+1}: '+' '.join(mp.nstr(v,3) for v in vs))
print('=== std T3 rows (weight 4, proved) ===')
for lab,row in std_small:
    terms=[(tuple(w),ab) for w,ab in row]
    vs=[abs(evalrel(terms,k)) for k in (1,3,5,7)]
    print(lab,' '.join(mp.nstr(v,3) for v in vs))
print('=== random control (must NOT vanish) ===')
random.seed(12345)
terms=[(tuple(fw3[j]),(0,random.randint(-1000,1000))) for j in range(32)]
for k in (1,3,5,7):
    print(f'k={k}: |v|={mp.nstr(abs(evalrel(terms,k)),4)}')
