"""Constant-retaining point symbol and a necessary classicality test.

The two primitive Lie words of length three are paired against the
antisymmetrized (3,2) deconcatenation. This kills shuffle products and Li5.
A zero answer is inconclusive: primitive zeta(3) is lost by the full symbol.
"""
import json
import sys
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'A3_W4_proof'))
from w4_symbol import *


def expression():
    out = {}
    for z, factors, coefficient in [
            (ZERO, (LT,LT,MT,MT), 7*R/128),
            (ZERO, (LT,LT,LT,MT), -5*R/96),
            (ZERO, (LT,LT,LT,LT), Q(F(1,24))-3*R/128),
            (-A, (LT,LT,LT,LT), Q(F(-1,12)))]:
        for (w,k), d in log_product(*factors).items():
            key = ((z,)+w,k)
            out[key] = out.get(key,ZERO)+coefficient*d
    return clean(out)


def symbol():
    out = {}
    for (w,k),c in expression().items():
        s = add(g_symbol(ONE,w,A),scale(g_symbol(ONE,w,B),-1))
        for _ in range(k):
            s = shuffle(s,log_symbol(Q(2)))
        out = add(out,scale(s,c))
    return out


def bracket(x,y):
    return add(concat(x,y),scale(concat(y,x),-1))


X,Y = {(0,):ONE},{(1,):ONE}
L2 = bracket(X,Y)
L3 = [bracket(X,L2),bracket(Y,L2)]


def obstruction(s):
    return [sum((c*d*(s.get(u+v,ZERO)-s.get(v+u,ZERO))
                 for u,c in p.items() for v,d in L2.items()),ZERO) for p in L3]


if __name__ == '__main__':
    ban_float()
    self_check()
    # The same test must kill every product split and classical Li5 atom.
    controls = 0
    for k in range(1,5):
        for u in product(range(2),repeat=k):
            for v in product(range(2),repeat=5-k):
                assert obstruction(shuffle({u:ONE},{v:ONE})) == [ZERO,ZERO]
                controls += 1
    for z in [B,-B,-B*B,R/2,Q(F(1,2)),ONE-R/2,2-R]:
        assert obstruction(li_symbol(5,z)) == [ZERO,ZERO]
        controls += 1
    s = symbol()
    obs = obstruction(s)
    data = dict(gpl_terms=len(expression()), product_and_Li5_controls=controls,
                alphabet=['log2','log(1+sqrt2)'],
                symbol={''.join(map(str,w)):[str(c.p),str(c.q)] for w,c in sorted(s.items())},
                delta32_symbol_pairing=[[str(c.p),str(c.q)] for c in obs])
    (HERE/'symbol.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2),flush=True)
    print('Nonzero is an obstruction; zero does NOT prove classicality.',flush=True)
