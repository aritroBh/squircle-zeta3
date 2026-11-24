"""Inherited slow expressions vs fast feeds; negative control perturbs a family input."""
from pathlib import Path
import types
import sympy as sp
from sympy.polys.fields import field
from sympy.polys.domains import QQ_I
import prove_adj910_fast as fast
src=Path('../TOW3_tower_depth6/prove_adj56.py').read_text()
src=src.replace('import tow3lib as T3','import tow5lib as T3').replace('T3.e(', 'T3.forced(').replace('T3.tau1(', 'T3.tau1f(').replace('T3.tau2(', 'T3.tau2f(')
slow=types.ModuleType('slow_expressions');exec(compile(src,'inherited_slow_expressions','exec'),slow.__dict__)
p=fast.p; K,P=field(p,QQ_I)
for dl,ss in [(9,3),(10,4)]:
 st=slow.feed_terms(dl,ss);ft=fast.feed_terms(dl,ss)
 assert len(st)==len(ft)
 for j,(a,b) in enumerate(zip(st,ft)):
  assert K.from_expr(a)==K.from_expr(b),(dl,ss,j)
 lhs=slow.eR(dl,ss,slow.d,p)+slow.eR(dl,ss+1,slow.d,p)
 residual=K.from_expr(lhs)-sum((K.from_expr(t) for t in st),K.zero)
 assert residual==K.zero
 print(f'Slow un-cached expressions vs fast ({dl},{ss}): lhs, all {len(st)} terms, and slow field residual PASS.',flush=True)
# Perturb an actual family input, rather than the final numerator.
original=fast.eR
def mutated(dl,s,dd,mm,conj=False):
 v=original(dl,s,dd,mm,conj)
 return v+1 if (dl,s)==(9,4) else v
fast.eR=mutated;fast._ehat.clear();fast._ed.clear()
try:
 fast.main(9,4)
except AssertionError:
 print('EXPECTED FAILURE: fast checker rejects F(9,4)+1 input mutation.',flush=True)
else:
 raise AssertionError('Mutation escaped checker')
fast.eR=original;fast._ehat.clear();fast._ed.clear();fast.main(9,4)
print('Restored input passes the same fast checker.',flush=True)
