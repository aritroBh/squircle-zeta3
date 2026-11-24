import sys, json, math, time, itertools
from fractions import Fraction as Ff
sys.path.insert(0,'.')
import w5blib as W, w5verify as V

gl = json.load(open('out/d13/genuine_lattice_cap.json'))
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
vnum=[int(x) for x in json.load(open('out/d13/gcd_vec_num.json'))['v_num']]
cur_shift=json.load(open('out/d13/gcd_vec_num.json'))['shift']
Ror=[omap(r) for r in L0o]; Obase=omap(vnum); m=len(Obase)
# NOTE: Obase corresponds to cur_shift already applied; wider box {-2..2}^4 short x {-1,0,1}^4 long AROUND vnum
t0=time.time()
best=(max(abs(x) for x in Obase), (0,)*8)
cnt=0
for s4 in itertools.product((-2,-1,0,1,2), repeat=4):
    O1=[Obase[t]+sum(s4[i]*Ror[i][t] for i in range(4)) for t in range(m)]
    for s8 in itertools.product((-1,0,1), repeat=4):
        O=[O1[t]+sum(s8[i]*Ror[4+i][t] for i in range(4)) for t in range(m)]
        sc=max(abs(x) for x in O); cnt+=1
        if sc<best[0]: best=(sc, s4+s8)
print(f'wider enum {cnt} ({time.time()-t0:.0f}s); best {best[0]} shift {best[1]} (base {max(abs(x) for x in Obase)})')
# coordinate descent from best
cur=list(best[1]); curs=best[0]; it=0; improved=True
while improved and it<300:
    improved=False; it+=1
    for i in range(8):
        for sgn in (1,-1):
            k=1
            while k<(1<<22):
                trial=cur[:]; trial[i]+=sgn*k
                O=[Obase[t]+sum(trial[j]*Ror[j][t] for j in range(8)) for t in range(m)]
                sc=max(abs(x) for x in O)
                if sc<curs: cur=trial; curs=sc; improved=True; k*=2
                else: break
print(f'descent2: score {curs} shift {cur} iters {it}')
if curs < max(abs(x) for x in Obase):
    vnum=[vnum[j]+sum(cur[i]*L0o[i][j] for i in range(8)) for j in range(n)]
    assert vnum[0]==g
e=expand(vnum)
coefs=list(e['li5'].values())+list(e['prod'].values())+[e['a_z5']]
mn=max(abs(c.numerator) for c in coefs); lm=1
for c in coefs: lm=math.lcm(lm,c.denominator)
print('FINAL v_num: nli5',len(e['li5']),'nprod',len(e['prod']),'maxnum',mn,'lcm',lm,'z5',e['a_z5'],'atom-max',max(abs(x) for x in vnum))
json.dump(dict(v_num=[str(x) for x in vnum], stats=dict(nli5=len(e['li5']),nprod=len(e['prod']),maxnum=str(mn),lcm=lm)),
          open('out/d13/gcd_vec_num.json','w'))
print('saved out/d13/gcd_vec_num.json')

# ---- sanity checks on BOTH candidate vectors (cheap): exact symbol, factorization, resid @300
from w5blib import QD
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),13), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
F=W.make_field(13, split=[3], inert=[2])
vatom=[int(x) for x in json.load(open('out/d13/gcd_vec_atom.json'))['v_atom']]
for tag,v in (('v_atom',vatom),('v_num',vnum)):
    ee=expand(v)
    t0=time.time(); sym=V.symbol_check5(cands, ee, names)
    nchk,bad=V.factorization_check(cands, ee, F)
    mp_resid=V.explicit_resid(13, {'D':13,'names':names,'letters':F['letters']}, cands, ee, 300)
    print(f'{tag}: symbol nonzero={len(sym)}; factor checked {nchk} bad {bad}; resid@300 {mp.nstr(mp_resid,3)} ({time.time()-t0:.0f}s)')
