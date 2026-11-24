"""Same submitted adjacent feeds; exact QQ(i)(p) arithmetic avoids huge DEN products."""
from pathlib import Path
import sys,time
import sympy as sp
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'TOW3_tower_depth6'))
import prove_adj56 as A
K=sp.QQ_I.frac_field(A.p)
for delta in (5,6):
    for s in range(delta//2+1):
        start=time.monotonic()
        lhs=A.eR(delta,s,A.d,A.p)+A.eR(delta,s+1,A.d,A.p)
        rem=K.from_sympy(lhs)
        for term in A.feed_terms(delta,s): rem-=K.from_sympy(term)
        assert rem==K.zero,(delta,s,rem)
        print(f'PASS adjacent ({delta},{s}) exact QQ(i)(p) remainder 0; {time.monotonic()-start:.2f}s',flush=True)
