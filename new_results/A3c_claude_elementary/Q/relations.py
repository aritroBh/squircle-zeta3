"""Theorem relations through weight 4; all matrix elimination is in python-flint.

H words as in cmzv.py. Extended DS: reg_shuffle(u*v-u shuffle v),
u in A1, v admissible. Regularized distribution uses Zhao (5.4).
No period-independence or completeness assumption is made.
"""
from cmzv import *
from math import gcd, lcm, factorial
from flint import nmod_mat, fmpq_mat, ctx
import argparse, pickle, signal
ctx.threads = 1
PRIME = 1000000009  # 1 mod 8; primality checked in main


@cache
def words(n,admissible=True):
    if not n: return ((),)
    return tuple(w for w in product(range(-1,8),repeat=n)
                 if w[-1] != -1 and (not admissible or w[0] != 0))


@cache
def to_series(w):
    out=[]; s=1; last=0
    for j in w:
        if j == -1: s+=1
        else:
            out.append((s,(j-last)%8)); last=j; s=1
    assert s == 1
    return tuple(out)


def from_series(s):
    out=[]; last=0
    for n,j in s:
        last=(last+j)%8
        out.extend([-1]*(n-1)+[last])
    return tuple(out)


@cache
def stuffle(u,v):
    if not u: return {v:1}
    if not v: return {u:1}
    out=Counter()
    for w,a in stuffle(u[1:],v).items(): out[(u[0],)+w]+=a
    for w,a in stuffle(u,v[1:]).items(): out[(v[0],)+w]+=a
    merged=(u[0][0]+v[0][0],(u[0][1]+v[0][1])%8)
    for w,a in stuffle(u[1:],v[1:]).items(): out[(merged,)+w]+=a
    return dict(out)


def ds(u,v):
    a={from_series(w):c for w,c in stuffle(to_series(u),to_series(v)).items()}
    return regularize(add((a,1),(sh(u,v),-1)))


def distribution(n,d):
    C={(j,):F(1) for j in range(8//d,8,8//d)}
    cp=[{():F(1)}]
    for k in range(1,n+1): cp.append(add((mul(cp[-1],C),F(1,k))))
    for w in product((-1,)+tuple(range(0,8,d)),repeat=n):
        if w[-1] == -1: continue
        lift=[((-1,),) if j == -1 else tuple((j//d+k*(8//d),) for k in range(d)) for j in w]
        left=Counter()
        for tup in product(*lift): left[tuple(x[0] for x in tup)]+=d**w.count(-1)
        right=reg(w)
        for k in range(1,n+1):
            if w[k-1] != 0: break
            right=add((right,1),(mul(cp[k],reg(w[k:])),1))
        yield regularize(add((dict(left),1),(right,-1)))


def primitive(row):
    row={w:F(c) for w,c in row.items() if c}
    if not row: return ()
    den=lcm(*(c.denominator for c in row.values()))
    vals=sorted((w,c.numerator*(den//c.denominator)) for w,c in row.items())
    g=gcd(*(c for w,c in vals))
    if vals[0][1] < 0: g=-g
    return tuple((w,c//g) for w,c in vals)


def make_relations(n,lower):
    rows=[]; seen=set(); counts=Counter()
    def emit(row,family):
        r=primitive(row)
        if not r: return
        counts[family+'_nonzero']+=1
        if r not in seen:
            seen.add(r); rows.append(r); counts[family+'_new']+=1
    if n==1:
        for j in range(1,8):
            emit(add((add(({(j,):1},1),({((-j)%8,):1},-1)),2),({(2,):1,(6,):-1},-(4-j))), 'weight1_branch')
    else:
        for a in range(1,n):
            b=n-a
            for u in words(a,False):
                u_adm=u[0]!=0
                if u_adm and a>b: continue
                for v in words(b):
                    if u_adm and a==b and u>v: continue
                    emit(ds(u,v),'double_shuffle')
    for d in (2,4,8):
        for row in distribution(n,d): emit(row,'distribution')
    for k,basis in lower.items():
        for row in basis:
            for w in words(n-k): emit(mul(dict(row),{w:1}),'lift')
    log(f'weight {n} relation families: {dict(counts)}; unique rows={len(rows)}')
    return rows,dict(counts)


def matrix(rows,cols,p=PRIME):
    ids={w:i for i,w in enumerate(cols)}
    mat=nmod_mat(len(rows),len(cols),p)
    for i,row in enumerate(rows):
        for w,c in row: mat[i,ids[w]]=c%p
    return mat


def pivots(rref,rank):
    out=[]; start=0
    for i in range(rank):
        while not rref[i,start]: start+=1
        out.append(start); start+=1
    return out


def small_calibration(n,rows,cols):
    start=time.monotonic()
    mat=matrix(rows,cols)
    rref,rank=mat.rref()
    # Pivot columns of transpose select actual theorem rows, preserving their
    # small rational coefficients for subsequent lifts and exact verification.
    trans,rank2=mat.transpose().rref()
    assert rank==rank2
    selected=[rows[j] for j in pivots(trans,rank)]
    qmat=fmpq_mat(len(selected),len(cols))
    ids={w:i for i,w in enumerate(cols)}
    for i,row in enumerate(selected):
        for w,c in row: qmat[i,ids[w]]=c
    qr,qrank=qmat.rref()
    assert qrank==rank
    # Verify ALL extra rows are rational consequences using the exact nullspace
    # of the selected row space. Dimension is small; multiplication is in FLINT.
    pp=pivots(qr,rank)
    free=[j for j in range(len(cols)) if j not in pp]
    null=fmpq_mat(len(cols),len(free))
    for k,j in enumerate(free):
        null[j,k]=1
        for i,p in enumerate(pp): null[p,k]=-qr[i,j]
    allq=fmpq_mat(len(rows),len(cols))
    for i,row in enumerate(rows):
        for w,c in row: allq[i,ids[w]]=c
    assert allq*null==fmpq_mat(len(rows),len(free))
    log(f'weight {n}: columns={len(cols)}, EXACT Q rank={rank}, quotient={len(free)}, DG benchmark={3**n}; seconds={time.monotonic()-start:.3f}')
    return selected,{'weight':n,'columns':len(cols),'relations':len(rows),'rank_Q':rank,'quotient':len(free),'DG':3**n},null


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--max-weight',type=int,default=3)
    ap.add_argument('--prime',type=int,default=PRIME)
    ap.add_argument('--kummer',action='store_true')
    args=ap.parse_args()
    from flint import fmpz
    assert fmpz(args.prime).is_prime() and args.prime%8==1
    log(f'Relation calibration started: prime={args.prime}, FLINT threads={ctx.threads}')
    lower={}; summary=[]
    for n in range(1,args.max_weight+1):
        cols=words(n); rows,counts=make_relations(n,lower)
        if n==4 and args.kummer:
            with (HERE/'kummer_aux.pickle').open('rb') as f: extra=pickle.load(f)[3]
            seen=set(rows); added=0
            for row in extra:
                r=primitive(row)
                if r and r not in seen: rows.append(r); seen.add(r); added+=1
            counts['kummer_new']=added
            log(f'Kummer auxiliary-elimination rows added={added}')
        if n<=3:
            selected,report,null=small_calibration(n,rows,cols)
            lower[n]=selected
            with (HERE/f'weight{n}.pickle').open('wb') as f: pickle.dump((cols,rows,selected),f)
        else:
            start=time.monotonic()
            with (HERE/'weight4.pickle').open('wb') as f: pickle.dump((cols,rows),f)
            log(f'weight 4 matrix allocation: {len(rows)} x {len(cols)}; raw data bytes={8*len(rows)*len(cols)}')
            mat=matrix(rows,cols,args.prime)
            log('weight 4 FLINT modular rref started')
            rr,rank=mat.rref(inplace=True)
            pp=pivots(rr,rank); free=[j for j in range(len(cols)) if j not in pp]
            log(f'weight 4: modular rank={rank}, quotient={len(free)}, DG benchmark={3**n}; seconds={time.monotonic()-start:.3f}')
            data=json.loads((HERE/'vectors.json').read_text())
            ids={w:i for i,w in enumerate(cols)}
            residues=[]
            for part in data['D']:
                v=[0]*len(cols)
                for w,c in part:
                    c=F(c); v[ids[tuple(w)]]=c.numerator*pow(c.denominator,-1,args.prime)%args.prime
                rem=[]
                for j in free:
                    rem.append((v[j]-sum(v[p]*int(rr[i,j]) for i,p in enumerate(pp)))%args.prime)
                log(f'weight 4 D residue nonzero coordinates={sum(bool(x) for x in rem)} of {len(free)}')
                residues.append(rem)
            report={'weight':n,'columns':len(cols),'relations':len(rows),'rank_mod_p':rank,'prime':args.prime,'quotient':len(free),'DG':3**n,'D_nonzero':[sum(bool(x) for x in r) for r in residues]}
            pmap={p:i for i,p in enumerate(pp)}
            (HERE/'modular_nullspace.json').write_text(json.dumps({'prime':args.prime,'free_words':[cols[j] for j in free], 'rows':[[(-int(rr[pmap[i],j]))%args.prime if i in pmap else int(i==j) for j in free] for i in range(len(cols))]}))
            (HERE/'modular_residue.json').write_text(json.dumps({'prime':args.prime,'free_words':[cols[j] for j in free],'residues':residues}))
        report['families']=counts; summary.append(report)
        (HERE/'calibration.json').write_text(json.dumps(summary,indent=2))
    log('Relation calibration completed.')

if __name__=='__main__':
    signal.alarm(3600)
    main()
