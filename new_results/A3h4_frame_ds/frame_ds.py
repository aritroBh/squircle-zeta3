"""Function double shuffle in real projective frames; single process, exact rows.

Li_s(z)=(-1)^depth G(0^(s1-1),1/z1,...;1).
The domain used here is real 0<|z_j|<1 for every input argument.
This implies absolute convergence of every stuffle and shuffle integral.
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
HERE = Path(__file__).resolve().parent
THEORY = HERE.parent
sys.path.insert(0, str(THEORY/'A3c_claude_elementary/Q'))
from cmzv import Q, F, POINTS, cross, dlog, forms_path, add, mul, sh, path
from relations import words, primitive, matrix
from itertools import permutations, product
from functools import cache
from collections import Counter
import pickle, json, time, signal, argparse
from flint import nmod_mat, ctx
ctx.threads = 1
START = time.monotonic()

def log(s):
    s = time.strftime('%Y-%m-%d %H:%M:%S')+' '+s
    print(s, flush=True)
    with (HERE/'progress.log').open('a') as f: f.write(s+'\n')

@cache
def compositions(n):
    if n == 0: return ((),)
    return tuple((k,)+t for k in range(1,n+1) for t in compositions(n-k))

@cache
def stuffle(u,v):
    if not u: return {v:1}
    if not v: return {u:1}
    out = Counter()
    for w,c in stuffle(u[1:],v).items(): out[(u[0],)+w] += c
    for w,c in stuffle(u,v[1:]).items(): out[(v[0],)+w] += c
    merged = (u[0][0]+v[0][0],u[0][1]*v[0][1])
    for w,c in stuffle(u[1:],v[1:]).items(): out[(merged,)+w] += c
    return dict(out)

def gword(u):
    out=[]; b=Q(1)
    for s,z in u:
        b *= z
        out.extend([Q(0)]*(s-1)+[1/b])
    return tuple(out)

def identity(factors):
    """Return exact G words for series stuffle minus integral shuffle."""
    st={():1}; sf={():1}; depth=sum(map(len,factors))
    for u in factors:
        st=add(*[(stuffle(v,u),c) for v,c in st.items()])
        sf=add(*[(sh(v,gword(u)),c) for v,c in sf.items()])
    return add(*[({gword(u):1},c*(-1)**len(u)) for u,c in st.items()],
               (sf,-(-1)**depth))

def frames():
    for p0,p1,pinf in permutations(POINTS,3):
        ts={cross(p,p0,p1,pinf):p for p in POINTS}
        assert len(ts)==len(POINTS) and ts[Q(0)]==p0 and ts[Q(1)]==p1 and ts[None]==pinf
        bs={1/t for t in ts if t is not None and abs(t)>1}
        yield (p0,p1,pinf),ts,bs

def sequences(bs,d):
    if not d: return ((),)
    return tuple(b for b in product(bs,repeat=d)
                 if all(abs(b[i])<abs(b[i-1]) for i in range(1,d)))

def factor(s,b):
    return tuple((k,x/(b[i-1] if i else Q(1))) for i,(k,x) in enumerate(zip(s,b)))

def candidates(n,bs,arity=2):
    """All ordered factors, modulo permutation of the factors.

    Every lattice vertex in a stuffle path occurs in some stuffle term;
    therefore all nonempty subset products of cumulative b's must be in bs.
    """
    for weights in compositions(n):
        if len(weights)!=arity: continue
        for ss in product(*(compositions(k) for k in weights)):
            for bb in product(*(sequences(bs,len(s)) for s in ss)):
                fs=tuple(factor(s,b) for s,b in zip(ss,bb))
                keys=[repr(f) for f in fs]
                if keys!=sorted(keys): continue
                if all(__import__('functools').reduce(lambda x,y:x*y,vertex,Q(1)) in bs
                       for vertex in product(*((Q(1),)+b for b in bb)) if any(x!=1 for x in vertex)):
                    yield fs

def convert(g,ts,rel):
    p0,p1,pinf=g
    return add(*[(forms_path([dlog(ts[a],pinf) for a in w],POINTS[p0],POINTS[p1]),c)
                 for w,c in rel.items()])

def dump(name,data):
    with (HERE/name).open('wb') as f: pickle.dump(data,f)

def load(path):
    with Path(path).open('rb') as f: return pickle.load(f)

def numerical_selfcheck():
    import mpmath as mp
    mp.mp.dps=65
    def mq(q): return mp.mpf(str(q.p))+mp.sqrt(2)*mp.mpf(str(q.q))
    def series(u,N=450):
        tail=[mp.mpf(1)]*(N+1)
        for s,z in reversed(u):
            z=mq(z); acc=mp.mpf(0); new=[mp.mpf(0)]*(N+1)
            for k in range(1,N+1):
                acc += z**k*tail[k-1]/k**s
                new[k]=acc
            tail=new
        return tail[N]
    def G(w,N=450):
        cs=[mp.mpf(1)]+[mp.mpf(0)]*N
        for a in reversed(w):
            new=[mp.mpf(0)]*(N+1)
            if a==0:
                assert not cs[0]
                for k in range(1,N+1): new[k]=cs[k]/k
            else:
                a=mq(a); v=mp.mpf(0)
                for k in range(N):
                    v=(v-cs[k])/a; new[k+1]=v/(k+1)
            cs=new
        return mp.fsum(cs)
    checks=[]
    for n in (2,3,4):
        for weights in compositions(n):
            if len(weights) not in (2,3,4): continue
            for ss in product(*(compositions(k) for k in weights)):
                fs=tuple(tuple((s,Q(F(1,4+i+j))) for j,s in enumerate(comp)) for i,comp in enumerate(ss))
                lhs=mp.fprod(series(f) for f in fs)
                st={():1}
                for f in fs: st=add(*[(stuffle(v,f),c) for v,c in st.items()])
                err=abs(lhs-mp.fsum(c*series(v) for v,c in st.items()))
                errg=abs(mp.fsum(c*G(w) for w,c in identity(fs).items()))
                checks.append(max(err,errg))
                assert max(err,errg)<mp.mpf('1e-40'),(fs,err,errg)
    log(f'direct mpmath series/G selfcheck instances={len(checks)} dps={mp.mp.dps} maximum residual={mp.nstr(max(checks),8)}')

def generate(n):
    rows={}; counts=Counter(); start=time.monotonic()
    for i,(g,ts,bs) in enumerate(frames()):
        counts['frames']+=1
        for arity in range(2,n+1):
            for fs in candidates(n,bs,arity):
                counts[f'candidates_{arity}']+=1
                rel=identity(fs)
                assert all(a in ts for w in rel for a in w)
                row=primitive(convert(g,ts,rel))
                if row: rows.setdefault(row,(g,fs))
                else: counts['zero']+=1
        if i%32==31:
            log(f'w{n} frames={i+1} unique_nonzero={len(rows)} seconds={time.monotonic()-start:.1f}')
        if time.monotonic()-start>2700: raise TimeoutError('generation step exceeded time limit')
    out={'rows':list(rows),'provenance':list(rows.values()),'counts':dict(counts)}
    dump(f'frame_w{n}.pickle',out)
    log(f'w{n} generation complete {dict(counts)} unique_nonzero={len(rows)} seconds={time.monotonic()-start:.2f}')
    return out

def validate(n,data):
    import mpmath as mp
    mp.mp.dps=100
    if n==3:
        raw=load(THEORY/'A3f6_level16/vals8_150digit.pickle')
        vals={w:mp.mpc(*v) for w,v in raw.items()}
    else:
        sys.path.insert(0,str(THEORY/'A3f6_level16'))
        import eval8 as ev
        from flint import acb,arb
        ctx.prec=400; ctx.cap=1000; ev.N=900
        ev.roots=[acb(0,arb.pi()*j/4).exp() for j in range(8)]
        ev.roots[0],ev.roots[2],ev.roots[4],ev.roots[6]=acb(1),acb(0,1),acb(-1),acb(0,-1)
        vals={}
        for w in set(w for row in data['rows'] for w,c in row):
            z=ev.h(w)
            vals[w]=mp.mpc(z.real.mid().str(110,radius=False,more=True),z.imag.mid().str(110,radius=False,more=True))
        log(f'w{n} numerical H words evaluated={len(vals)} precision_bits={ctx.prec} cap={ctx.cap} terms={ev.N}')
        dump(f'values_w{n}.pickle',{w:(str(z.real),str(z.imag)) for w,z in vals.items()})
    accepted=[]; provenance=[]; rejected=[]; worst=mp.mpf(0)
    for row,prov in zip(data['rows'],data['provenance']):
        residual=abs(mp.fsum(mp.mpf(c)*vals[w] for w,c in row))
        scale=max(1,sum(abs(c) for w,c in row))
        if residual/scale<mp.mpf('1e-65' if n==3 else '1e-45'):
            accepted.append(row); provenance.append(prov); worst=max(worst,residual)
        else: rejected.append((row,prov,str(residual)))
    data.update(rows=accepted,provenance=provenance,rejected=rejected,worst=str(worst))
    dump(f'validated_w{n}.pickle',data)
    log(f'w{n} numerical accepted={len(accepted)} rejected={len(rejected)} max_accepted_absolute_residual={mp.nstr(worst,8)}')
    return data

def rank2(data):
    cols,base,selected=load(THEORY/'A3c_claude_elementary/Q/weight2.pickle')
    out=[]
    for p in (1000000009,1000000033):
        b=matrix(base,cols,p).rank(); r=matrix(base+data['rows'],cols,p).rank()
        log(f'w2 p={p} columns={len(cols)} base_rank={b} plus_frames_rank={r} quotient={len(cols)-r} target_rank={len(cols)-3**2}')
        out.append(dict(p=p,columns=len(cols),base=b,rank=r,quotient=len(cols)-r))
    (HERE/'ranks_w2.json').write_text(json.dumps(out,indent=2))
    return out

PRIMES=((1000000009,291087696),(1000000033,12641440))

@cache
def fmod(c,p):
    c=F(c); return c.numerator*pow(c.denominator,-1,p)%p

def base3():
    b=load(THEORY/'A3f3_level8_compute/weight3_with_lifts.pickle')
    # The file's producer step_51.py stores (q,p) = p+q*sqrt(2).
    k=load(THEORY/'A3f3_level8_compute/k3_newHrows.pickle')['krows']
    return [tuple(c) for c in b['cols']], [tuple((tuple(w),(0,c)) for w,c in row) for row in b['rows']]+[tuple(a.items()) for a in k]

def kmatrix(rows,cols,p,r):
    ids={w:i for i,w in enumerate(cols)}; M=nmod_mat(len(rows),len(cols),p)
    for i,row in enumerate(rows):
        for w,(q,c) in row: M[i,ids[w]]=(fmod(q,p)*r+fmod(c,p))%p
    return M

def pivotcols(rr,rank):
    out=[]; j=0
    for i in range(rank):
        while not rr[i,j]: j+=1
        out.append(j); j+=1
    return out

def rank3(data):
    import mpmath as mp
    mp.mp.dps=90
    vals={w:mp.mpc(*v) for w,v in load(THEORY/'A3f6_level16/vals8_150digit.pickle').items()}
    cols,base=base3(); extras=[tuple((w,(0,c)) for w,c in row) for row in data['rows']]
    k=load(THEORY/'A3f3_level8_compute/k3_newHrows.pickle')['krows']
    for swapped in (False,True):
        err=max(abs(mp.fsum((mp.mpf(str(y))+mp.sqrt(2)*mp.mpf(str(x)) if not swapped else mp.mpf(str(x))+mp.sqrt(2)*mp.mpf(str(y)))*vals[w] for w,(x,y) in a.items())) for a in k)
        log(f'K3 coefficient orientation swapped={swapped} max numerical residual={mp.nstr(err,8)}')
    five=load(THEORY/'A3f4_glanois_basis/five_relations.pickle')
    targets=[tuple((tuple(w),(0,c)) for w,c in zip(five['freewords'],rel) if c) for rel in five['relations']]
    results=[]
    for p,r in PRIMES:
        b=kmatrix(base,cols,p,r).rank(); M=kmatrix(base+extras,cols,p,r)
        rr,rank=M.rref(); pp=pivotcols(rr,rank); free=[j for j in range(len(cols)) if j not in pp]
        nf=kmatrix(targets,cols,p,r)
        members=[]
        for k in range(len(targets)):
            res=[(int(nf[k,j])-sum(int(nf[k,h])*int(rr[i,j]) for i,h in enumerate(pp)))%p for j in free]
            members.append(not any(res))
        log(f'w3 p={p} sqrt2={r} columns={len(cols)} base_rank={b} plus_frames_rank={rank} quotient={len(cols)-rank} target_rank={len(cols)-3**3} R_members={members}')
        results.append(dict(p=p,r=r,base=b,rank=rank,quotient=len(cols)-rank,R_members=members))
    (HERE/'ranks_w3.json').write_text(json.dumps(results,indent=2))
    return results

def null4(saved,p):
    cols=words(4); pp=saved['pp']; free=saved['free']
    assert saved['rank']==len(pp) and sorted(pp+free)==list(range(len(cols)))
    ids={j:k for k,j in enumerate(free)}
    N=nmod_mat(len(cols),len(free),p)
    for j,k in ids.items(): N[j,k]=1
    for h,row in zip(pp,saved['rowsF']):
        for j,c in row: N[h,ids[j]]=-c%p
    return N

def reduce4(rows,N,p,r):
    return kmatrix(rows,words(4),p,r)*N

def rank4(data):
    start=time.monotonic(); cols=words(4)
    extras=[tuple((w,(0,c)) for w,c in row) for row in data['rows']]
    # Include all proved new lower-weight identities, not just selected directions.
    for n in (2,3):
        for row in load(HERE/f'validated_w{n}.pickle')['rows']:
            for w in words(4-n):
                extras.append(tuple((t,(0,c)) for t,c in mul(dict(row),{w:1}).items()))
    log(f'w4 own-weight frame rows={len(data["rows"])} plus all lower lifts total_extra_rows={len(extras)}')
    targets=load(THEORY/'A3h_weight4_target/w4_systems.pkl')
    systems=load(THEORY/'A3h_weight4_target/w4_systems_T.pkl')
    reports=[]
    for system in ('T3','T3noR'):
        for p,r in PRIMES:
            saved=load(THEORY/f'A3h_weight4_target/rref_{system}_{p}.pkl')
            assert saved['p']==p and saved['r']==r
            N=null4(saved,p)
            # Verify every source row annihilates the supplied normal-form map.
            for i in range(0,len(systems[system]),128):
                batch=[row for lab,row in systems[system][i:i+128]]
                reduced=reduce4(batch,N,p,r)
                assert reduced==nmod_mat(len(batch),N.ncols(),p),(system,p,i)
            log(f'w4 {system} p={p} all_source_rows={len(systems[system])} checked against supplied RREF nullspace; cached_base_rank={saved["rank"]}')
            nr=reduce4(extras,N,p,r); rk=nr.rank()
            own=reduce4(extras[:len(data['rows'])],N,p,r).rank()
            pp4=targets['P4']; D={}
            for w,(q,c) in targets['Drat'].items():
                assert q==0; D[w]=(0,c)
            for w,(q,c) in targets['Dsq'].items():
                assert q==0
                a,b=D.get(w,(0,0)); D[w]=(a+c,b)
            v=reduce4([tuple(pp4.items()),tuple(D.items())],N,p,r)
            withp=nmod_mat(nr.tolist()+[v.tolist()[0]],p).rank()
            withd=nmod_mat(nr.tolist()+v.tolist(),p).rank()
            report=dict(system=system,p=p,r=r,base_rank=saved['rank'],own_added_rank=own,
                        added_rank=rk,rank=saved['rank']+rk,quotient=len(cols)-saved['rank']-rk,
                        extra_plus_P4_rank=withp,extra_plus_P4_D_rank=withd,D_member=withd==withp)
            reports.append(report)
            log('w4 '+json.dumps(report))
            if time.monotonic()-start>2700: raise TimeoutError('weight-four analysis exceeded time limit')
    (HERE/'ranks_w4.json').write_text(json.dumps(reports,indent=2))
    return reports

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('action',choices=['selfcheck','generate','validate','rank']); ap.add_argument('weight',type=int,nargs='?',default=2)
    args=ap.parse_args(); signal.alarm(2700)
    if args.action=='selfcheck': numerical_selfcheck()
    elif args.action=='generate': generate(args.weight)
    elif args.action=='validate': validate(args.weight,load(HERE/f'frame_w{args.weight}.pickle'))
    elif args.action=='rank':
        {2:rank2,3:rank3,4:rank4}[args.weight](load(HERE/f'validated_w{args.weight}.pickle'))
