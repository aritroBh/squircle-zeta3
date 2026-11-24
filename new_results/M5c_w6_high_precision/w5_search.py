"""Exact point-symbol elimination followed by a smaller, exact-Gram LLL.

No generic-function symbol is used. MPL convention: Li_(s,t)(x,y) =
sum_(n>m>0) x^n y^m/(n^s m^t). Two selected Li_(3,2) atoms span the
visible obstruction; additional atoms test its full-symbol kernel.
"""
import itertools
import json
import math
import signal
import sys
import time
from fractions import Fraction as F
from pathlib import Path

import flint
from flint import acb, arb, ctx
import mpmath as mp
import w5_symbol as W
from w5_precision import log

HERE = Path(__file__).resolve().parent
ZERO, ONE, Q, R, A, B = W.ZERO,W.ONE,W.Q,W.R,W.A,W.B


def sf(s, c):
    return W.scale(s,Q(c))


def generate(extra, m3_basis=False):
    # Each item is (name, symbol, numerical recipe).
    groups = {
        1:[('lam',W.log_symbol(Q(2)),('log',Q(2))),
           ('L',W.log_symbol(A),('log',A))],
        2:[('pi2',{},('pi2',)),('Li2(b)',W.li_symbol(2,B),('li',2,B))],
        3:[('z3',{},('zeta',3))]+[(f'Li3({z})',W.li_symbol(3,z),('li',3,z))
                                for z in [B,-B,R/2]],
        4:[(f'Li4({z})',W.li_symbol(4,z),('li',4,z)) for z in
           [Q(F(1,2)),B,-B,R/2,-B*B,2-R,ONE-R/2]]}
    if m3_basis:
        # Actual greedy survivors printed in m3_pslq_1150.out, not M4a's fixed-list typo.
        groups[4]=[(f'Li4({z})',W.li_symbol(4,z),('li',4,z)) for z in
                   [Q(F(1,2),F(1,4)),-2+2*R,R/2,-R/2,2-R,B,-B]]
        log('Using the actual seven greedy Li4 survivors from m3_pslq_1150.out')
    out = []
    for partition in [(1,1,1,1,1),(2,1,1,1),(2,2,1),(3,1,1),(3,2),(4,1)]:
        choices = [list(itertools.combinations_with_replacement(groups[k],partition.count(k)))
                   for k in sorted(set(partition))]
        for blocks in itertools.product(*choices):
            atoms = sum(blocks,())
            s = {():ONE}
            for _,v,_ in atoms:
                s = W.shuffle(s,v)
            out.append(('*'.join(a[0] for a in atoms),s,('product',tuple(a[2] for a in atoms))))
    log(f'Weight-five product count={len(out)}')
    candidates = set()
    for j in range(-6,7):
        for m in range(-4,5):
            for sign in [-1,1]:
                z = Q(sign)*R**j*A**m
                if not -ONE <= z < ONE or z == ZERO:
                    continue
                try:
                    W.logvec(ONE-z)
                except ValueError:
                    continue
                candidates.add(z)
    for z in sorted(candidates,key=lambda z:(z.p,z.q)):
        out.append((f'Li5({z})',W.li_symbol(5,z),('li',5,z)))
    out.append(('z5',{},('zeta',5)))
    extras = [(3,2,x,y) for x,y in [(B,ONE),(R/2,ONE),(-B,ONE),(B,-ONE),(R/2,-ONE),
              (ONE-R/2,ONE),(2-R,ONE),(-B*B,ONE)]]
    extras += [(4,1,x,ONE) for x in [B,-B,R/2,ONE-R/2,2-R]]
    for s,t,x,y in extras[:extra]:
        word = (ZERO,)*(s-1)+(ONE/x,)+(ZERO,)*(t-1)+(ONE/(x*y),)
        out.append((f'Li{s}{t}({x},{y})',W.g_symbol(ZERO,word,ONE),('mpl',s,t,x,y)))
    log(f'Li5 S-unit candidates={len(candidates)}; depth-two atoms={min(extra,len(extras))}; total atoms={len(out)}')
    return out


def eliminate(basis,target):
    # Rational symbols: the target has two rational prefactor parts.
    words = list(itertools.product(range(2),repeat=5))
    assert all(c.q == 0 for _,s,_ in basis for c in s.values())
    rows = [[flint.fmpq(str(s.get(w,ZERO).p)) for _,s,_ in basis]
            +[flint.fmpq(str(target.get(w,ZERO).p)),flint.fmpq(str(target.get(w,ZERO).q))]
            for w in words]
    mat,rank = flint.fmpq_mat(rows).rref()
    pivots = [next(j for j in range(mat.ncols()) if mat[i,j]) for i in range(rank)]
    assert all(j<len(basis) for j in pivots), 'extra atoms fail to span target symbol'
    free = [j for j in range(len(basis)) if j not in pivots]
    correction = [(p,Q(F(str(mat[i,len(basis)])),F(str(mat[i,len(basis)+1]))))
                  for i,p in enumerate(pivots)]
    kernels = [(j,[(p,F(str(mat[i,j]))) for i,p in enumerate(pivots) if mat[i,j]]) for j in free]
    # Independent reconstruction, including both algebraic prefactor parts.
    s = target
    for j,c in correction:
        s = W.add(s,W.scale(basis[j][1],-c))
    assert not s
    for j,cs in kernels:
        s = basis[j][1]
        for p,c in cs:
            s = W.add(s,sf(basis[p][1],-c))
        assert not s
    log(f'Exact symbol rank={rank}; remaining field-valued kernel atoms={len(kernels)}')
    for j,c in correction:
        if basis[j][0].startswith('Li32') and c != ZERO:
            log(f'Exact obstruction correction: {c} * {basis[j][0]}')
    return correction,kernels


def numeric(z):
    return arb(z.p.numerator)/z.p.denominator+arb(2).sqrt()*arb(z.q.numerator)/z.q.denominator


def mpl(s,t,x,y,digits):
    # ponytail: convergent |x|<1, |y|<=1 series only; use analytic continuation for other arguments.
    xx,yy = numeric(x),numeric(y)
    r = abs(float(xx))
    assert 0<r<1 and abs(yy)<=1
    N = math.ceil((digits+25)*math.log(10)/-math.log(r))+10
    total,H,power,yp = arb(0),arb(0),arb(1),arb(1)
    for n in range(1,N+1):
        power *= xx
        total += power*H/(n**s)
        yp *= yy
        H += yp/(n**t)
    rr = abs(xx)
    tail = rr**(N+1)*((N+1)/(1-rr)+rr/((1-rr)*(1-rr)))
    assert tail < arb(10)**(-digits)
    return total+arb(0,tail.abs_upper())


def evaluate(recipe,digits):
    kind,*a = recipe
    if kind == 'log': return numeric(a[0]).log()
    if kind == 'pi2': return arb.pi()*arb.pi()
    if kind == 'zeta': return arb(a[0]).zeta()
    if kind == 'li':
        z = acb(numeric(a[1])).polylog(a[0])
        assert z.imag.contains(0)
        return z.real
    if kind == 'mpl': return mpl(*a,digits)
    assert kind == 'product'
    ans = arb(1)
    for atom in a[0]: ans *= evaluate(atom,digits)
    return ans


def lll(vals,digits,label):
    t0 = time.monotonic()
    n = len(vals)
    C = mp.mpf(10)**digits
    rows = [[int(i==j) for j in range(n)]+[int(mp.nint(C*v))] for i,v in enumerate(vals)]
    log(f'LLL start {label}: rows={n}, digits={digits}, gram=exact')
    signal.alarm(2700)
    reduced = flint.fmpz_mat(rows).lll(gram='exact')
    signal.alarm(0)
    relations = []
    for i in range(n):
        vec = [int(reduced[i,j]) for j in range(n)]
        mx = max(map(abs,vec))
        residual = mp.fsum(c*v for c,v in zip(vec,vals))
        # Require unused decimal digits as well as a short coefficient vector.
        # A short row at lattice precision alone can be a random approximation.
        if abs(residual)<mp.mpf(10)**(-digits-30) and mx<mp.mpf(10)**(mp.mpf(digits)/n):
            relations.append(dict(vec=vec,maxcoef=mx,residual=mp.nstr(residual,12)))
    log(f'LLL done {label}: seconds={time.monotonic()-t0:.3f}; small relations={len(relations)}; first maxcoef={max(abs(int(reduced[0,j])) for j in range(n))}')
    return relations


if __name__ == '__main__':
    digits = int(sys.argv[1]) if len(sys.argv)>1 else 1200
    mode = sys.argv[2] if len(sys.argv)>2 else 'coaction'
    extra = 0 if mode.startswith('coaction') else int(mode)
    # Precision AFTER all imports; do not reuse roots created at lower precision.
    ctx.prec = math.ceil((digits+120)*math.log2(10))
    mp.mp.dps = digits+100
    test = mpl(1,1,B,ONE,digits+20)-numeric(ONE-B).log()**2/2
    assert test.contains(0)
    log(f'MPL series control Li11(b,1)=log(1-b)^2/2: {test.str(8)}')
    plant = [mp.mpf(1),mp.sqrt(2),1+mp.sqrt(2)]
    planted = lll(plant,digits,'planted relation')
    assert any(r['vec'] in [[1,1,-1],[-1,-1,1]] for r in planted)
    basis = generate(extra,mode=='coaction_m3')
    target_symbol=W.symbol()
    fixed=[]
    if mode.startswith('coaction'):
        certificate=json.loads((HERE/'cobracket.json').read_text())
        for term in certificate['correction']:
            x,y,c=(Q(*map(F,term[k])) for k in ['x','y','coefficient'])
            word=(ZERO,ZERO,ONE/x,ZERO,ONE/(x*y))
            target_symbol=W.add(target_symbol,W.scale(W.g_symbol(ZERO,word,ONE),-c))
            fixed.append((c,('mpl',3,2,x,y)))
    correction,kernels = eliminate(basis,target_symbol)
    numbers = [evaluate(recipe,digits+60) for _,_,recipe in basis]
    maxrad = max(x.rad() for x in numbers)
    assert maxrad<arb(10)**(-digits-30)
    log(f'Atom evaluation maximum radius={maxrad.str(8)}')
    vals = [mp.mpf(v.mid().str(digits+90,radius=False)) for v in numbers]
    r2 = mp.sqrt(2)
    def qnum(c): return mp.mpf(c.p.numerator)/c.p.denominator+r2*mp.mpf(c.q.numerator)/c.q.denominator
    saved = json.loads((HERE/'W5_4400.json').read_text())
    if (HERE/'W5_5400.json').exists():
        saved=json.loads((HERE/'W5_5400.json').read_text())
    assert digits+80 < len(saved['value']), 'compute a higher-precision W5 before increasing lattice precision'
    fixed_value=mp.fsum(qnum(c)*mp.mpf(evaluate(recipe,digits+60).mid().str(digits+90,radius=False))
                       for c,recipe in fixed)
    target = mp.mpf(saved['value'])-fixed_value-mp.fsum(qnum(c)*vals[j] for j,c in correction)
    kv = [vals[j]-mp.fsum(mp.mpf(c.numerator)/c.denominator*vals[p] for p,c in cs) for j,cs in kernels]
    # Remove pre-existing basis dependencies, preserving explicit exact rational vectors.
    kept = list(range(2*len(kv)))
    values = lambda: [target]+[(r2 if j%2 else 1)*kv[j//2] for j in kept]
    relations = lll(values(),digits,'symbol-kernel enlarged basis')
    dependencies=[]
    while any(r['vec'][0]==0 for r in relations):
        deps=[r for r in relations if r['vec'][0]==0]
        dependencies.append(dict(kept=kept.copy(),relations=deps))
        mat,rank=flint.fmpq_mat([list(reversed(r['vec'][1:])) for r in deps]).rref()
        drop={len(kept)-1-next(j for j in range(mat.ncols()) if mat[i,j]) for i in range(rank)}
        kept=[j for i,j in enumerate(kept) if i not in drop]
        log(f'Removed {len(drop)} independent numerical basis dependencies; retained scalar atoms={len(kept)}')
        relations=lll(values(),digits,'independent symbol-kernel enlarged basis')
    data = dict(digits=digits,extra_atoms=extra,mode=mode,names=[n for n,_,_ in basis],
                correction=[(j,str(c.p),str(c.q)) for j,c in correction],
                kernels=[(j,[(p,str(c)) for p,c in cs]) for j,cs in kernels],
                kept=kept,relations=relations,dependencies=dependencies)
    (HERE/f'search_{digits}_{mode}.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2),flush=True)
