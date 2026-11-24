import sys, json, time
from fractions import Fraction as Ff
sys.path.insert(0,'.')
import w5blib as W, w5pipe
from w5blib import QD
# S = {29,2,5,7}: add split prime 7
F=W.make_field(29, split=[5,7], inert=[2])
print("letters:", F['names'], "d=", len(F['names']), flush=True)
for B in (4,5):
    t0=time.time()
    c = W.enum_exceptional_fast(F, B=B, time_limit=1500)
    from collections import Counter
    hc = Counter(x['h'] for x in c)
    print(f"B={B}: {len(c)} exceptional, heights {dict(sorted(hc.items()))} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(dict(D=29, B=B, split=[5,7], inert=[2], names=F['names'],
                   disp=[W.fmt_qd(g) for g in F['letters']], S=[29,2,5,7],
                   cands=[dict(x=[str(x['x'].a), str(x['x'].b)], vx=x['vx'], v1m=x['v1m'], h=x['h'],
                               xf=float(x['x'].emb(40))) for x in c]),
              open(f'out/d29/candidates_S7_B{B}.json','w'))
    print(f"saved out/d29/candidates_S7_B{B}.json", flush=True)
