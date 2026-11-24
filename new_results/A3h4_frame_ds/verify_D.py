"""Standalone exact Fraction verifier. Rebuilds identities and W4-R from sources.

No numerical values, modular membership conclusions, Kummer rows or R1--R5
enter the identity check. Existing lower-basis files select rows only: every
selected row is independently regenerated from standard theorem identities.
"""
from frame_ds import *
import cmzv as cm
import relations as standard

def rebuild_standard():
    standard.log=log  # All logs stay in this work directory.
    selected={}; maps={}
    for n in (1,2,3):
        generated,_=standard.make_relations(n,selected)
        oldcols,oldrows,basis=load(THEORY/f'A3c_claude_elementary/Q/weight{n}.pickle')
        assert tuple(oldcols)==words(n)
        generated=set(generated)
        assert all(row in generated for row in basis),('lower basis not regenerated',n)
        selected[n]=basis
        maps[n]=basis
    generated,_=standard.make_relations(4,{})
    return generated,maps

def rebuild_target():
    lhs=[{},{}]
    for (w,k),coefficient in cm.source.w4_g_expression().items():
        forms=[dlog(p) for p in w]
        term=add((forms_path(forms,POINTS[Q(1)],POINTS[cm.A]),1),
                 (forms_path(forms,POINTS[Q(1)],POINTS[cm.B]),-1))
        for _ in range(k): term=mul(term,cm.atom(('log',Q(2))))
        for j,c in enumerate((coefficient.p,coefficient.q)):
            lhs[j]=add((lhs[j],1),(term,c))
    rhs=[{},{}]
    for coefficient,atoms in cm.source.conjecture_terms():
        term={():F(1)}
        for a in atoms: term=mul(term,cm.atom(a))
        for j,c in enumerate((coefficient.p,coefficient.q)):
            rhs[j]=add((rhs[j],1),(term,c))
    return [add((lhs[j],1),(rhs[j],-1)) for j in range(2)]

def main():
    cert=json.loads((HERE/'D_certificate.json').read_text())
    std,bases=rebuild_standard()
    ipi={(2,):F(2),(6,):F(-2)}
    P4=mul(mul(ipi,ipi),mul(ipi,ipi))
    rebuilt={}; counts=Counter()
    for label in cert['labels']:
        key=repr(label)
        if label[0]=='T0':
            lab=label[1]; family=lab[0]
            if family=='std': row=std[lab[1]]
            else:
                assert family in ('lift1','lift2old','lift3old'),lab
                n={'lift1':1,'lift2old':2,'lift3old':3}[family]
                row=primitive(mul(dict(bases[n][lab[1]]),{tuple(lab[2]):1}))
            counts[family]+=1
        elif label[0]=='frame':
            _,n,i,w=label; data=load(HERE/f'validated_w{n}.pickle')
            g,fs=data['provenance'][i]
            p0,p1,pinf=g
            ts={cross(p,p0,p1,pinf):p for p in POINTS}
            bs={1/t for t in ts if t is not None and abs(t)>1}
            # Regenerate and validate the complete stuffle grid in exact Q(sqrt2).
            assert fs in tuple(candidates(n,bs,len(fs)))
            raw=primitive(convert(g,ts,identity(fs)))
            assert raw==data['rows'][i],('frame export mismatch',n,i)
            row=primitive(mul(dict(raw),{tuple(w):1}))
            counts[f'frame_w{n}']+=1
        else:
            assert label==['P4']; row=tuple(P4.items()); counts['P4']+=1
        rebuilt[key]=dict(row)
    targets=rebuild_target()
    c=[]
    for j,(target,coeffs) in enumerate(zip(targets,cert['coefficients'])):
        residual=Counter({w:-a for w,a in target.items()}); pi=F(0)
        for label,a in zip(cert['labels'],coeffs):
            a=F(a)
            if label==['P4']: pi=a
            for w,b in rebuilt[repr(label)].items(): residual[w]+=a*b
        nonzero={w:a for w,a in residual.items() if a}
        log(f'EXACT rebuilt D component={j} support={len(target)} residual_nonzero={len(nonzero)} c{j}={pi}')
        assert not nonzero,('exact identity failed',j,next(iter(nonzero.items()),None))
        c.append(pi)
    log(f'EXACT certificate source families={dict(counts)} total_used_rows={len(cert["labels"])}')
    log(f'VERDICT: D=({c[0]}+sqrt(2)*{c[1]})*pi^4; Conjecture 4 '+('PROVED' if c==[0,0] else 'NOT established'))

if __name__=='__main__': main()
