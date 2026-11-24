
import sys, os, pickle, time
from fractions import Fraction as F
A3H='<project>/research/08_theory/A3h_weight4_target'
sys.path.insert(0, os.path.join(A3H,'Q'))
from relations import words
from flint import nmod_mat
b12 = pickle.load(open(os.path.join(A3H,'basis_w12.pkl'),'rb'))
ws = pickle.load(open(os.path.join(A3H,'w4_systems.pkl'),'rb'))
P4, Drat, Dsq = ws['P4'], ws['Drat'], ws['Dsq']
Ts = pickle.load(open(os.path.join(A3H,'w4_systems_T.pkl'),'rb'))
cols = words(4); ids={w:i for i,w in enumerate(cols)}; NC=len(cols)
_modp_cache={}
def modp(c,p):
    key=(c,p)
    v=_modp_cache.get(key)
    if v is None:
        c=F(c); v=(c.numerator%p)*pow(c.denominator%p,-1,p)%p; _modp_cache[key]=v
    return v
def vec_mod(vdict,p,r):
    v=[0]*NC
    for w,(q,pp) in vdict.items():
        v[ids[w]]=(modp(q,p)*r+modp(pp,p))%p
    return v
def analyze(p,r,tag):
    rows=Ts['T3']
    t0=time.monotonic()
    M=nmod_mat(len(rows),NC,p)
    for i,(lab,row) in enumerate(rows):
        for w,(q,pp) in row:
            M[i,ids[w]]=((q%p)*r+(pp%p))%p
    t1=time.monotonic()
    RR,rank=M.rref(inplace=True)
    t2=time.monotonic()
    pset=set(); pp_=[]; start=0
    for i in range(rank):
        j=start
        while not RR[i,j]: j+=1
        pp_.append(j); pset.add(j); start=j+1
    free=[j for j in range(NC) if j not in pset]
    rowsF=[[(j,int(RR[i,j])) for j in free if RR[i,j]] for i in range(rank)]
    colpos={j:k for k,j in enumerate(free)}
    colview=[[] for _ in free]
    for i,rowf in enumerate(rowsF):
        for j,v in rowf: colview[colpos[j]].append((i,v))
    def nf(v):
        vp=[v[j]%p for j in pp_]
        out=[]
        for j,entries in zip(free,colview):
            s=v[j]%p
            for i,val in entries:
                vi=vp[i]
                if vi: s=(s-vi*val)%p
            out.append(s)
        return out
    res={}
    for name,vd in [('P4',P4),('Drat',Drat),('Dsq',Dsq)]:
        res[name]=nf(vec_mod(vd,p,r))
    print(f'  p={p} r={r}: rank={rank} quotient={len(free)} build={t1-t0:.1f}s rref={t2-t1:.1f}s nf-nz='+','.join(f'{n}={sum(bool(x) for x in res[n])}' for n in res),flush=True)
    pickle.dump({'sys':'T3','p':p,'r':r,'rank':rank,'quotient':len(free),'pp':pp_,'free':free,'rowsF':rowsF,'nf':res},
                open(f'<project>/research/08_theory/A3h2_w4_missing/rref_T3_{p}_{tag}.pkl','wb'))
if __name__=='__main__':
    p=int(sys.argv[1]); tag=sys.argv[2]
    import sympy as sp
    roots=sp.sqrt_mod(2,p,all_roots=True)
    r = roots[0] if tag=='a' else roots[1]
    print('prime',p,'root',r,tag)
    analyze(p,r,tag)
