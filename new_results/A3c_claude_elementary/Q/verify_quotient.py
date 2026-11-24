"""Reconstruct a modular nullspace; verify it exactly with FLINT, then print formal normal forms."""
from relations import *
from flint import fmpq_mat, fmpq
from math import isqrt


@cache
def reconstruct(a,m):
    bound=isqrt(m//2)
    r0,r1,t0,t1=m,a,0,1
    while abs(r1)>bound:
        q=r0//r1
        r0,r1,t0,t1=r1,r0-q*r1,t1,t0-q*t1
    if t1<0: r1,t1=-r1,-t1
    if not t1 or t1>bound or gcd(r1,t1)!=1 or (a*t1-r1)%m:
        raise ValueError(f'No bounded rational reconstruction: a={a}, modulus={m}')
    return fmpq(r1,t1)


def main():
    data=json.loads((HERE/'modular_nullspace.json').read_text())
    with (HERE/'weight4.pickle').open('rb') as f: cols,rows=pickle.load(f)
    p=data['prime']; raw=data['rows']; d=len(data['free_words'])
    log(f'Exact nullspace reconstruction started: rows={len(raw)}, columns={d}, modulus={p}')
    N=fmpq_mat(len(cols),d)
    try:
        for i,row in enumerate(raw):
            for j,x in enumerate(row): N[i,j]=reconstruct(x,p)
    except ValueError as e:
        log(str(e)); return False
    log(f'Nullspace entries reconstructed; distinct residues={reconstruct.cache_info().currsize}')
    ids={w:i for i,w in enumerate(cols)}
    # ponytail: small row batches cap memory; all elimination/multiplication is FLINT.
    start=time.monotonic()
    for offset in range(0,len(rows),128):
        batch=rows[offset:offset+128]
        M=fmpq_mat(len(batch),len(cols))
        for i,row in enumerate(batch):
            for w,c in row: M[i,ids[w]]=c
        residual=M*N
        if residual!=fmpq_mat(len(batch),d):
            log(f'Exact nullspace verification FAILED in batch offset={offset}; need more CRT primes')
            return False
        if offset and offset%2048==0: log(f'Exact relation rows verified={offset}; elapsed seconds={time.monotonic()-start:.3f}')
    log(f'ALL exact rational relation rows verified={len(rows)}; elapsed seconds={time.monotonic()-start:.3f}')
    # Identity on free words fixes normalization and proves nullspace independence.
    for j,w in enumerate(data['free_words']):
        for k in range(d): assert N[ids[tuple(w)],k]==int(j==k)
    vec=json.loads((HERE/'vectors.json').read_text())['D']
    pi4=mul(atom(('pi2',)),atom(('pi2',)))
    targets=fmpq_mat(3,len(cols))
    for i,part in enumerate(vec):
        for w,c in part: targets[i,ids[tuple(w)]]=c
    for w,c in pi4.items(): targets[2,ids[w]]=str(c)
    nf=targets*N
    exactrank=len(cols)-d
    log(f'CERTIFIED Q rank={exactrank}; exact quotient dimension={d}; DG benchmark={3**4}')
    _,targetrank=nf.rref()
    log(f'Formal span rank of (D0,D1,pi^4) normal forms={targetrank}')
    out={'prime_or_CRT_modulus':p,'rank_Q':exactrank,'quotient':d,
         'free_words':data['free_words'],
         'normal_forms':[[str(nf[i,j]) for j in range(d)] for i in range(3)],
         'target_span_rank':targetrank}
    (HERE/'exact_normal_form.json').write_text(json.dumps(out,indent=2))
    text=[]
    for i,label in enumerate(('D0','D1','pi^4')):
        terms=[f'{nf[i,j]} * H{tuple(w)}' for j,w in enumerate(data['free_words']) if nf[i,j]]
        log(f'{label}: exact formal nonzero normal-form coordinates={len(terms)}')
        text.append(label+' =\n'+'\n'.join(terms))
    (HERE/'normal_form.txt').write_text('\n\n'.join(text)+'\n')
    print('\n\n'.join(text),flush=True)
    return True

if __name__=='__main__':
    signal.alarm(3600)
    if not main(): sys.exit(1)
