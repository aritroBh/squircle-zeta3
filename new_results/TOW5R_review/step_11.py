"""step_11 (TOW5R): fast lcm-method vs SLOW plain sp.together method at depths 9, 10.
Recompute lhs - sum(terms) for adjacent (9,3),(9,4),(10,5) with plain together+cancel
(no DEN/lcm machinery): numerator must be exactly 0. Validates the fast method at the new depths.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/tow5_copy')
import sympy as sp
import importlib.util
spec = importlib.util.spec_from_file_location('paf', os.path.dirname(os.path.abspath(__file__)) + '/tow5_copy/prove_adj910_fast.py')
paf = importlib.util.module_from_spec(spec)
import types
# prevent main() from running: spec loader executes module top-level only (main guarded)
spec.loader.exec_module(paf)
p = paf.p
d = 2*p-1
ok_all = True
for (dl, s) in [(9,3),(9,4),(10,5)]:
    t0 = time.time()
    lhs = sp.cancel(paf.eR(dl, s, d, p) + paf.eR(dl, s+1, d, p))
    terms = paf.feed_terms(dl, s)
    tot = sp.Integer(0)
    for tm in terms:
        tot += tm
    rem = sp.together(lhs - tot)
    num = sp.expand(sp.numer(rem))
    ok = (num == 0)
    ok_all &= ok
    print(f'SLOW ({dl},{s}): numerator == 0 : {ok} ({round(time.time()-t0,1)}s, {len(terms)} terms)', flush=True)
print('STEP_11 VERDICT:', 'PASS' if ok_all else 'FAIL')
assert ok_all
