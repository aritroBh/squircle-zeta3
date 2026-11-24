"""Expand and independently check the enlarged-basis relation; write its formula."""
import json
import math
from collections import Counter
from fractions import Fraction as F
from functools import reduce
from pathlib import Path

import mpmath as mp
import sympy as sp
from flint import ctx
import w5_symbol as W
import w5_cobracket as C
from w5_search import generate,evaluate,numeric
from w5_precision import log

HERE=Path(__file__).resolve().parent


def expand(data):
    rs=[r for r in data['relations'] if r['vec'][0]]
    assert len(rs)==1
    relation=rs[0]
    vec=relation['vec']
    out=[W.ZERO]*len(data['names'])
    for j,p,q in data['correction']:
        out[j]+=W.Q(F(p),F(q))
    for v,k in zip(vec[1:],data['kept']):
        j,cs=data['kernels'][k//2]
        cf=W.Q(F(-v,vec[0]))*(W.R if k%2 else W.ONE)
        out[j]+=cf
        for p,c in cs: out[p]-=cf*W.Q(F(c))
    return out,relation


def texatom(recipe):
    kind,*a=recipe
    if kind=='log': return r'\lambda' if a[0]==W.Q(2) else 'L'
    if kind=='pi2': return r'\pi^2'
    if kind=='zeta': return rf'\zeta({a[0]})'
    if kind=='li': return rf'\operatorname{{Li}}_{a[0]}\!\left({sp.latex(a[1].sympy())}\right)'
    if kind=='mpl':
        s,t,x,y=a
        return rf'\operatorname{{Li}}_{{{s},{t}}}\!\left({sp.latex(x.sympy())},{sp.latex(y.sympy())}\right)'
    assert kind=='product'
    parts=[]
    for atom,n in Counter(a[0]).items():
        tex=texatom(atom)
        if n>1:
            if tex in ['L',r'\lambda']: tex+=rf'^{n}'
            elif atom[0]=='pi2': tex=rf'\pi^{{{2*n}}}'
            else: tex=rf'\left[{tex}\right]^{{{n}}}'
        parts.append(tex)
    return '\\,'.join(parts)


if __name__=='__main__':
    low=json.loads((HERE/'search_1200_13.json').read_text())
    high=json.loads((HERE/'search_1500_13.json').read_text())
    lc,lr=expand(low)
    hc,hr=expand(high)
    assert low['names']==high['names'] and lc==hc, 'relation not stable'
    log(f'Expanded relation identical at lattice digits={low["digits"]},{high["digits"]}; raw maxcoef={lr["maxcoef"]},{hr["maxcoef"]}')
    basis=generate(13)
    assert [n for n,_,_ in basis]==high['names']
    den=math.lcm(*(f.denominator for q in hc for f in [q.p,q.q]))
    ints=[den]+[int(-f*den) for q in hc for f in [q.p,q.q]]
    gcd=reduce(math.gcd,ints)
    ints=[i//gcd for i in ints]
    log(f'Original-basis primitive integer relation: W5 coefficient={ints[0]}, maxabs={max(map(abs,ints))}, nonzero field terms={sum(q!=W.ZERO for q in hc)}')

    # Check the complete point symbol, without numerical arithmetic.
    symbol={}
    for cf,(_,s,_) in zip(hc,basis):
        symbol=W.add(symbol,W.scale(s,cf))
    assert symbol==W.symbol()
    log('Expanded R5 full point symbol equals W5 exactly: PASS')

    # Independent numerical method for every depth-two atom actually used.
    mp.mp.dps=90
    ctx.prec=500
    for cf,(name,_,recipe) in zip(hc,basis):
        if cf==W.ZERO or recipe[0]!='mpl': continue
        _,s,t,x,y=recipe
        def qnum(q): return mp.mpf(q.p.numerator)/q.p.denominator+mp.sqrt(2)*mp.mpf(q.q.numerator)/q.q.denominator
        xx,yy=qnum(x),qnum(y)
        def integrand(u):
            if not u: return mp.mpf(0)
            li=-mp.log1p(-xx*yy*u) if t==1 else mp.polylog(t,xx*yy*u)
            return xx*(-mp.log(u))**(s-1)*li/(1-xx*u)/math.factorial(s-1)
        quad=mp.quad(integrand,[0,mp.mpf(1)/4,mp.mpf(1)/2,1])
        series=mp.mpf(evaluate(recipe,120).mid().str(140,radius=False))
        diff=abs(quad-series)
        log(f'Independent MPL quadrature {name}: absdiff={mp.nstr(diff,12)}')
        assert diff<mp.mpf('1e-80')

    # A fresh high-precision evaluation of the expanded, original-basis formula.
    ctx.prec=6200
    mp.mp.dps=1800
    value=mp.mpf(0)
    for cf,(_,_,recipe) in zip(hc,basis):
        if cf==W.ZERO: continue
        term=evaluate(recipe,1740)*numeric(cf)
        value+=mp.mpf(term.mid().str(1790,radius=False))
    saved=json.loads((HERE/'W5_5400.json').read_text())
    residual=value-mp.mpf(saved['value'])
    log(f'Expanded R5 minus saved W5: {mp.nstr(residual,16)}')
    assert abs(residual)<mp.mpf('1e-1500')

    # Full point cobracket, including its primitive zeta(3) directions.
    C.ban_float()
    right=C.P()
    for cf,(name,_,recipe) in zip(hc,basis):
        if cf==W.ZERO or recipe[0]!='mpl': continue
        _,s,t,x,y=recipe
        word=(W.ZERO,)*(s-1)+(W.ONE/x,)+(W.ZERO,)*(t-1)+(W.ONE/(x*y),)
        right+=C.delta_word(W.ZERO,word,W.ONE)*cf
    delta=C.reduce(right-C.target())
    assert not delta.d
    log(f'Expanded R5 full point (3,2) Lie cobracket minus W5: {delta}')

    output=dict(lattice_digits=[low['digits'],high['digits']],
                stable=True,raw_maxcoef=[lr['maxcoef'],hr['maxcoef']],
                integer_W5_coefficient=ints[0],integer_maxcoef=max(map(abs,ints)),
                residual=mp.nstr(residual,20),agreement_digits=1500,
                terms={name:[str(cf.p),str(cf.q)] for cf,(name,_,_) in zip(hc,basis) if cf!=W.ZERO})
    (HERE/'closed_form.json').write_text(json.dumps(output,indent=2)+'\n')
    lines=['# W5: candidate depth-two closed form','',
           'Let $b=\\sqrt2-1$, $L=\\log(1+\\sqrt2)$, and $\\lambda=\\log2$.',
           'The convention is $\\operatorname{Li}_{s,t}(x,y)=\\sum_{n>m>0}x^ny^m/(n^sm^t)$.','',
           'This is a numerical identity, stable at the two lattice precisions; its full point symbol and (3,2) Lie cobracket agree exactly. It is not an exact identity certificate.','',
           r'\[',r'\begin{aligned}',r'W_5={}&']
    count=0
    for cf,(_,_,recipe) in zip(hc,basis):
        if cf==W.ZERO: continue
        negative=cf<W.ZERO
        cf=-cf if negative else cf
        coeff='' if cf==W.ONE else sp.latex(cf.sympy())
        if cf.p and cf.q: coeff=rf'\left({coeff}\right)'
        sign='-' if negative else ('+' if count else '')
        lines.append(sign+coeff+texatom(recipe)+r'\\ &')
        count+=1
    lines[-1]=lines[-1].removesuffix(r'\\ &')
    lines += [r'\end{aligned}',r'\]','',
              f'Largest coefficient in the primitive integer relation: `{output["integer_maxcoef"]}`; coefficient of W5: `{ints[0]}`.',
              f'Fresh expanded-form residual against the saved quadrature string: `{output["residual"]}`.','',
              'All coefficients of this full expression are numerical relation-search results. The coaction correction in result.md is exact; it is a different choice of representatives modulo the kernel.']
    (HERE/'closed_form.md').write_text('\n'.join(lines)+'\n')
    log('Closed form and verification metadata written to closed_form.md and closed_form.json')
