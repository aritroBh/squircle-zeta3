"""Sparse modular certificates, CRT and rational reconstruction; no guessed rows."""
from frame_ds import *
from math import isqrt

def certificate_rows():
    src=load(THEORY/'A3h_weight4_target/w4_systems_T.pkl')['T0']
    rows=[]
    for lab,row in src:
        assert all(q==0 for w,(q,c) in row)
        rows.append((('T0',lab),tuple((w,c) for w,(q,c) in row)))
    for n in (4,2,3):
        data=load(HERE/f'validated_w{n}.pickle')
        for i,row in enumerate(data['rows']):
            for w in words(4-n):
                lifted=primitive(mul(dict(row),{w:1}))
                rows.append((('frame',n,i,w),lifted))
    v=load(THEORY/'A3h_weight4_target/w4_systems.pkl')
    rows.append((('P4',),tuple((w,c) for w,(q,c) in v['P4'].items())))
    targets=[tuple((w,c) for w,(q,c) in v[name].items()) for name in ('Drat','Dsq')]
    return rows,targets

def eliminate(rows,targets,p,selected=None):
    start=time.monotonic(); ids={w:i for i,w in enumerate(words(4))}
    order=selected if selected is not None else sorted(range(len(rows)),key=lambda i:(len(rows[i][1]),i))
    basis={}; recipes=[]; chosen=[]
    for count,src in enumerate(order):
        d={ids[w]:fmod(c,p) for w,c in rows[src][1] if fmod(c,p)}; trace=[]
        while d:
            j=min(d)
            if j not in basis: break
            b,k=basis[j]; a=d[j]; trace.append((k,a))
            for h,c in b.items():
                v=(d.get(h,0)-a*c)%p
                if v: d[h]=v
                else: d.pop(h,None)
        if d:
            j=min(d); inv=pow(d[j],-1,p); d={h:c*inv%p for h,c in d.items()}
            k=len(chosen); basis[j]=(d,k); chosen.append(src); recipes.append((inv,trace))
        elif selected is not None: raise ValueError(('bad prime lost selected pivot',p,src))
        if count%500==499:
            log(f'certificate p={p} processed={count+1} pivots={len(chosen)} elapsed={time.monotonic()-start:.1f}s')
        if time.monotonic()-start>2700: raise TimeoutError('certificate modular solve exceeded time limit')
    answers=[]
    for target in targets:
        d={ids[w]:fmod(c,p) for w,c in target if fmod(c,p)}; acc=[0]*len(chosen)
        while d:
            j=min(d)
            if j not in basis: raise ValueError(('target not member',p,j))
            b,k=basis[j]; a=d[j]; acc[k]=(acc[k]+a)%p
            for h,c in b.items():
                v=(d.get(h,0)-a*c)%p
                if v: d[h]=v
                else: d.pop(h,None)
        out=[0]*len(chosen)
        for k in reversed(range(len(chosen))):
            inv,trace=recipes[k]; a=acc[k]*inv%p; out[k]=a
            if a:
                for h,c in trace: acc[h]=(acc[h]-a*c)%p
        answers.append(out)
    # Independent modular check against original rows (not the echelon form).
    for target,answer in zip(targets,answers):
        check=Counter({w:-fmod(c,p) for w,c in target})
        for src,a in zip(chosen,answer):
            if a:
                for w,c in rows[src][1]: check[w]=(check[w]+a*fmod(c,p))%p
        assert all(c%p==0 for c in check.values())
    log(f'certificate p={p} original-row modular identity verified selected={len(chosen)} seconds={time.monotonic()-start:.2f}')
    return chosen,answers

def reconstruct(a,m):
    bound=isqrt(m//2); r0,r1=m,a; t0,t1=0,1
    while abs(r1)>bound:
        q=r0//r1; r0,r1=r1,r0-q*r1; t0,t1=t1,t0-q*t1
    if t1==0 or abs(t1)>bound: return None
    c=F(r1,t1)
    return c if (c.numerator-a*c.denominator)%m==0 else None

def exact_check(rows,targets,selected,coeffs):
    for target,cs in zip(targets,coeffs):
        check=Counter({w:-F(c) for w,c in target})
        for i,a in zip(selected,cs):
            if a:
                for w,c in rows[i][1]: check[w]+=a*c
        if any(check.values()): return False
    return True

def main():
    from sympy import nextprime
    rows,targets=certificate_rows(); selected=None; residues=None; modulus=1
    p=PRIMES[0][0]
    while True:
        cachefile=HERE/f'certificate_mod_{p}.pickle'
        if cachefile.exists():
            d=load(cachefile); chosen,answers=d['selected'],d['answers']
            assert selected is None or selected==chosen
        else:
            chosen,answers=eliminate(rows,targets,p,selected)
            dump(cachefile.name,dict(selected=chosen,answers=answers))
        selected=chosen
        if residues is None: residues=answers
        else:
            inv=pow(modulus,-1,p)
            residues=[[a+modulus*((b-a)*inv%p) for a,b in zip(old,new)] for old,new in zip(residues,answers)]
        modulus*=p
        coeffs=[[reconstruct(a,modulus) for a in part] for part in residues]
        failed=sum(c is None for part in coeffs for c in part)
        log(f'CRT modulus_bits={modulus.bit_length()} reconstruction_failed={failed}')
        if not failed and exact_check(rows,targets,selected,coeffs):
            used=[k for k in range(len(selected)) if any(part[k] for part in coeffs)]
            cert=dict(labels=[rows[selected[k]][0] for k in used],
                      coefficients=[[str(part[k]) for k in used] for part in coeffs],
                      modulus=str(modulus),primes_used=[])
            (HERE/'D_certificate.json').write_text(json.dumps(cert,indent=2))
            log(f'EXACT D certificate PASS used_rows={len(used)} families={dict(Counter(str(x[0]) for x in cert["labels"]))}')
            for j,part in enumerate(coeffs):
                c=next((part[k] for k,i in enumerate(selected) if rows[i][0]==('P4',)),F(0))
                log(f'EXACT c{j}={c} maximum_coefficient_height={max(max(abs(a.numerator),a.denominator) for a in part)}')
            break
        if not failed: log('reconstructed identity failed exact check; adding another prime')
        p=int(nextprime(p))
        while p%8!=1: p=int(nextprime(p))

if __name__=='__main__': main()
