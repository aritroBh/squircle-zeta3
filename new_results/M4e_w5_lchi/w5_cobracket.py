"""Point-level (3,2) Lie cobracket, including primitive trilogarithms.

Project weight two onto Li2(b) modulo products and pi^2, using its injective
antisymmetric symbol. Modulo products, shuffle gives
G_c(x,y,z;e) = - int_c^e log|v-x| log|v-z| dv/(v-y).
The supplied, generalized weight-three primitives evaluate this integral.
"""
import json
import sys
from functools import cache
from pathlib import Path

import w5_symbol as W
import flint
from w5_precision import log
from w4_right_coproduct import (P,PE,ZERO,ONE,Q,R,A,B,F,logp,pairprimitive,
                                li3_reduce,ban_float,modpi)

HERE = Path(__file__).resolve().parent


def indecomp(p):
    return P({k:c for k,c in p.d.items() if len(k)==1 and
              (k==('z3',) or isinstance(k[0],tuple) and k[0][0]=='Li3')})


@cache
def g2(c,w,e):
    s = W.g_symbol(c,w,e)
    return (s.get((0,1),ZERO)-s.get((1,0),ZERO))*2


@cache
def g3(c,w,e):
    if c==e: return P()
    x,y,z = w
    p = indecomp(-(pairprimitive(x,z,y,e)-pairprimitive(x,z,y,c)))
    # Check both primitive Lie pairings, which kill the discarded products.
    s = W.g_symbol(c,w,e)
    t = {}
    for (atom,),cf in p.d.items():
        if atom=='z3': continue
        assert atom[0]=='Li3'
        t = W.add(t,W.scale(W.li_symbol(3,atom[1]),cf))
    for lie in W.L3:
        assert sum((v*(s.get(u,ZERO)-t.get(u,ZERO)) for u,v in lie.items()),ZERO)==ZERO,(c,w,e,p)
    return p


def delta_word(c,w,e):
    assert len(w)==5
    points=(e,)+w+(c,)
    ans=P()
    for k in [2,3]:
        for i in range(6-k):
            gap=w[i:i+k]
            retained=w[:i]+w[i+k:]
            ci,ei=points[i+k+1],points[i]
            if k==2:
                ans += g3(c,retained,e)*g2(ci,gap,ei)
            else:
                ans -= g3(ci,gap,ei)*g2(c,retained,e)
    return ans


def reduce(p):
    return indecomp(li3_reduce(p))


def target():
    ans=P()
    for (w,k),c in W.expression().items():
        if k==0:
            ans += (delta_word(ONE,w,A)-delta_word(ONE,w,B))*c
    return reduce(ans)


def depth2(x,y):
    return reduce(delta_word(ZERO,(ZERO,ZERO,ONE/x,ZERO,ONE/(x*y)),ONE))


if __name__=='__main__':
    ban_float()
    for z in [B,-B,R/2,ONE-R/2,2-R]:
        red=reduce(g3(ZERO,(ZERO,ZERO,ONE/z),ONE)+PE.LI3(z))
        assert not red.d, (z,red)
    log('Weight-three projection versus classical Li3, including primitive constants: PASS')
    # Classical Li5 control, plus a nontrivial product control.
    assert not reduce(delta_word(ZERO,(ZERO,ZERO,ZERO,ZERO,ONE/B),ONE)).d
    product=P()
    for w in W.shuffles((ZERO,ONE/B),(ZERO,ZERO,ONE/B)):
        product += delta_word(ZERO,w,ONE)
    assert not reduce(product).d
    log('Full point cobracket classical Li5 and shuffle-product controls: PASS')
    lhs=target()
    log(f'Full point delta32(W5), tensor Li2(b): {lhs}')
    out={'W5':str(lhs),'depth_two':{}}
    selected=[]
    for x,y in [(B,ONE),(R/2,ONE),(-B,ONE),(B,-ONE),(R/2,-ONE),
                (ONE-R/2,ONE),(2-R,ONE),(-B*B,ONE)]:
        red=depth2(x,y)
        name=f'Li32({x},{y})'
        out['depth_two'][name]=str(red)
        log(f'Full point delta32({name}), tensor Li2(b): {red}')
        if len(selected)<4: selected.append((x,y,red))
    atoms=[('z3',),(('Li3',B),),(('Li3',-B),),(('Li3',R/2),)]
    matrix=flint.fmpq_mat([[str(red.d.get(at,ZERO).p) for _,_,red in selected] for at in atoms])
    assert all(c.q==0 for _,_,red in selected for c in red.d.values())
    rhs=flint.fmpq_mat([[str(lhs.d.get(at,ZERO).p),str(lhs.d.get(at,ZERO).q)] for at in atoms])
    sol=matrix.inv()*rhs
    correction=[]
    residual=lhs
    for i,(x,y,red) in enumerate(selected):
        c=Q(F(str(sol[i,0])),F(str(sol[i,1])))
        residual-=red*c
        correction.append(dict(x=[str(x.p),str(x.q)],y=[str(y.p),str(y.q)],
                               coefficient=[str(c.p),str(c.q)]))
        log(f'Exact full cobracket correction: {c} * Li32({x},{y})')
    assert not residual.d
    out['correction']=correction
    log(f'Full cobracket of W5 minus the four-atom correction: {residual}')
    (HERE/'cobracket.json').write_text(json.dumps(out,indent=2)+'\n')
