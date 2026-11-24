"""Real-domain Kummer Li4 instances whose every argument is a configuration cross ratio."""
from cmzv import *
from relations import primitive


def args(x,y):
    xi=1-x; eta=1-y
    return [(-x*x*y*eta/xi,1),(-y*y*x*xi/eta,1),
            (x*x*y/(eta*eta*xi),1),(y*y*x/(xi*xi*eta),1),
            (x*y,-6),(x*y/(xi*eta),-6),(-x*y/xi,-6),(-x*y/eta,-6),
            (x*eta,-3),(y*xi,-3),(-x*eta/xi,-3),(-y*xi/eta,-3),
            (-x/(xi*eta),-3),(-y/(xi*eta),-3),(x/eta,-3),(y/xi,-3),
            (x,6),(y,6),(-x/xi,6),(-y/eta,6)]


def instances():
    config={}
    for p0,p1,pinf,e in permutations(POINTS,4):
        z=cross(e,p0,p1,pinf)
        if z is not None and z<=Q(1):
            try: config[z]=find_cross(z)
            except ValueError: pass
    positive=sorted(z for z in config if Q(0)<z<Q(1))
    log(f'Kummer search: real cross ratios <=1={len(config)}; positive candidates={len(positive)}')
    good=[]; excluded=Counter()
    for i,x in enumerate(positive):
        for y in positive[i:]:
            if x+y>Q(1): continue
            aa=args(x,y)
            missing=[z for z,c in aa if z not in config and z not in (Q(0),Q(1))]
            if missing:
                excluded['outside_configuration']+=1
                continue
            row={}
            for z,c in aa:
                if z==Q(0): continue
                li={(-1,-1,-1,0):F(1)} if z==Q(1) else atom(('li',4,z))
                row=add((row,1),(li,c))
            logs=mul(atom(('li',1,x)),atom(('li',1,x)))
            logs=mul(logs,mul(atom(('li',1,y)),atom(('li',1,y))))
            row=add((row,1),(logs,F(-3,2)))
            good.append((x,y,row))
    log(f'Kummer admissible instances={len(good)}; nonzero formal rows={sum(bool(r) for x,y,r in good)}; exclusions={dict(excluded)}')
    return good


if __name__=='__main__':
    import pickle
    good=instances()
    with (HERE/'kummer4.pickle').open('wb') as f: pickle.dump(good,f)
    for x,y,row in good:
        log(f'Kummer x={x}, y={y}; exact formal support={len(row)}')


def auxiliary_closure():
    """Eliminate nonconfiguration Li4 atoms from real D0 Kummer instances.

    This is a finite search, not a completeness assertion. All auxiliary
    constants cancel by exact rational row operations before conversion to H.
    """
    from flint import fmpq_mat
    from math import comb
    config={z for p0,p1,pinf,e in permutations(POINTS,4)
            if (z:=cross(e,p0,p1,pinf)) is not None and z<=Q(1)}|{Q(1)}
    candidates=set()
    for j in range(-6,7):
        for m in range(-6,7):
            x=SQ2**j*A**m
            if not Q(0)<x<Q(1): continue
            try: source.logvec(1-x)
            except ValueError: continue
            candidates.add(x)
    candidates=sorted(candidates)
    log(f'Auxiliary Kummer search: positive S-unit candidates={len(candidates)}; exponents bounded by 6')
    def logpoly(z):
        a,b=source.logvec(z)
        return {(1,0,0):a,(0,1,0):b}
    def pmul(a,b):
        out=Counter()
        for k,c in a.items():
            for l,d in b.items(): out[tuple(x+y for x,y in zip(k,l))]+=c*d
        return {k:c for k,c in out.items() if c}
    def power(a,n):
        out={(0,0,0):F(1)}
        for _ in range(n): out=pmul(out,a)
        return out
    def li(z):
        if z==Q(0): return {}
        if z==Q(1): return {(0,0,2):F(1,90)}
        if z==Q(-1): return {(0,0,2):F(-7,720)}
        if z < Q(-1):
            logs=logpoly(-z)
            out={1/z:F(-1),(0,0,2):F(-7,360)}
            out=add((out,1),(power(logs,4),F(-1,24)),
                    (pmul(power(logs,2),{(0,0,1):F(1)}),F(-1,12)))
            return out
        return {z:F(1)}
    relations=[]; xy=[]
    for i,x in enumerate(candidates):
        for y in candidates[i:]:
            if x+y>Q(1): continue
            row={}
            for z,c in args(x,y): row=add((row,1),(li(z),c))
            row=add((row,1),(pmul(power(logpoly(1-x),2),power(logpoly(1-y),2)),F(-3,2)))
            relations.append(row); xy.append((x,y))
    atoms={z for row in relations for z in row if isinstance(z,Q)}
    added=0
    for z in sorted(atoms):
        if Q(0)<z<=Q(1) and z*z in atoms:
            relations.append(add((li(z*z),1),(li(z),-8),(li(-z),-8))); added+=1
    keys={z for row in relations for z in row}
    aux=sorted(z for z in keys if isinstance(z,Q) and z not in config)
    good=sorted(z for z in keys if isinstance(z,Q) and z in config)
    logs=sorted(z for z in keys if isinstance(z,tuple))
    cols=aux+good+logs; ids={z:i for i,z in enumerate(cols)}
    mat=fmpq_mat(len(relations),len(cols))
    for i,row in enumerate(relations):
        for z,c in row.items(): mat[i,ids[z]]=str(c)
    rr,rank=mat.rref()
    outputs=[]
    for i in range(rank):
        if any(rr[i,j] for j in range(len(aux))): continue
        hrow={}
        for j,z in enumerate(good,len(aux)):
            if not rr[i,j]: continue
            conv=atom(('li',4,z)) if z!=Q(1) else {(-1,-1,-1,0):F(1)}
            hrow=add((hrow,1),(conv,F(str(rr[i,j]))))
        for j,k in enumerate(logs,len(aux)+len(good)):
            if not rr[i,j]: continue
            conv={():F(1)}
            for a,n in zip((('log',Q(2)),('log',A),('pi2',)),k):
                for _ in range(n): conv=mul(conv,atom(a))
            hrow=add((hrow,1),(conv,F(str(rr[i,j]))))
        outputs.append(hrow)
    log(f'Auxiliary Kummer: D0 instances={len(xy)}, duplication rows={added}, auxiliary atoms={len(aux)}, configuration atoms={len(good)}, log monomials={len(logs)}, Q rank={rank}, eliminated rows={len(outputs)}, nonzero H rows={sum(bool(r) for r in outputs)}')
    import pickle
    with (HERE/'kummer_aux.pickle').open('wb') as f: pickle.dump((xy,relations,cols,outputs),f)
    return outputs
