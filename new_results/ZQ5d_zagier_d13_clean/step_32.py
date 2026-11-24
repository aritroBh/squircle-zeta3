import sys, json, math, time, itertools
from fractions import Fraction as Ff
sys.path.insert(0,'.')
import flint
import w5blib as W, w5verify as V

gl = json.load(open('out/d13/genuine_lattice_cap.json'))
G = [[int(x) for x in row] for row in gl['G']]
atom_names = gl['atom_names']; n=len(atom_names)
J = json.load(open('out/d13/candidates_cap.json'))
Kk = json.load(open('out/d13/kernel5_cap.json'))
K=Kk['K']; comps=Kk['comps']; names=J['names']
keepc=[0,1,2,3]+list(range(23,51)); g=36771482407; L=120
def canon_key(parts): return tuple(sorted(parts))
saved=json.load(open('out/d13/gcd_L0_ometric.json'))
L0o=[[int(x) for x in r] for r in saved['L0o']]
prodkeys=[tuple(k.split('*')) for k in saved['prodkeys']]; PK={k:i for i,k in enumerate(prodkeys)}
def expand(v): return V.expand_a(13, J['cands'], K, comps, keepc, atom_names, v, names)
def omap(v):
    e=expand(v)
    li5=[0]*len(J['cands'])
    for i,c in e['li5'].items(): li5[i]=int(c)*L
    pr=[0]*len(prodkeys)
    for k,c in e['prod'].items(): pr[PK[canon_key(k)]]+=int(c*L)
    return li5+pr+[int(e['a_z5'])*L]
vnum=[int(x) for x in json.load(open('out/d13/gcd_vec_num_babai.json'))['v_num']]
assert vnum[0]==g
t0=time.time()
Ror=[omap(r) for r in L0o]   # reduced O-images
Obase=omap(vnum)
m=len(Obase)
print(f'O images rebuilt ({time.time()-t0:.0f}s), m={m}, base max|O|={max(abs(x) for x in Obase)}')

def score(O): return max(abs(x) for x in O)
t0=time.time()
best=(score(Obase), (0,)*8)
cnt=0
for shift in itertools.product((-1,0,1), repeat=8):
    O=[Obase[t]+sum(shift[i]*Ror[i][t] for i in range(8)) for t in range(m)]
    s=score(O); cnt+=1
    if s<best[0]: best=(s,shift)
print(f'enum {cnt} ({time.time()-t0:.0f}s); best score {best[0]} shift {best[1]}')
# coordinate descent with growing steps from the best
cur=list(best[1]); curs=best[0]
improved=True; it=0
while improved and it<200:
    improved=False; it+=1
    for i in range(8):
        # try steps +-1,+-2,+-4,... +-2^20
        for sgn in (1,-1):
            k=1
            while k<(1<<21):
                trial=cur[:]; trial[i]+=sgn*k
                O=[Obase[t]+sum(trial[j]*Ror[j][t] for j in range(8)) for t in range(m)]
                s=score(O)
                if s<curs: cur=trial; curs=s; improved=True; k*=2
                else: break
print(f'coord descent: final score {curs} at {cur} (iters {it})')
vfin=[vnum[j]+sum(cur[i]*L0o[i][j] for i in range(8)) for j in range(n)]
assert vfin[0]==g
e=expand(vfin)
def stats(e):
    coefs=list(e['li5'].values())+list(e['prod'].values())+[e['a_z5']]
    mn=max(abs(c.numerator) for c in coefs); lm=1
    for c in coefs: lm=math.lcm(lm,c.denominator)
    return dict(nli5=len(e['li5']), nprod=len(e['prod']), maxnum=mn, lcm=lm, z5=str(e['a_z5']))
st=stats(e)
print('v_fin:', st, 'atom-max:', max(abs(x) for x in vfin), 'maxnum/L check:', curs//L, curs%L==0)
json.dump(dict(v_num=[str(x) for x in vfin], shift=cur, stats=st), open('out/d13/gcd_vec_num.json','w'))
print('saved out/d13/gcd_vec_num.json')
