import pickle, sys
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
base2 = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
rw = pickle.load(open(base2+"rewrite_exact.pkl","rb"))
need = set()
for r in rw.values():
    need.update(r.keys())
print("words appearing in rewrite:", len(need))
# closure under contiguous subword + mrev (letter map x -> 1-x on QS2: 0->1,1->0,p->q,q->p)
import oct5d2_engine as E
ZERO, ONE, P, Q = E.ZERO, E.ONE, E.P_LET, E.Q_LET
def mrev(v): return tuple(ONE-x for x in reversed(v))
closure = set(need)
changed = True
while changed:
    changed = False
    for w in list(closure):
        for i in range(len(w)+1):
            for u in (w[:i], w[i:]):
                if u not in closure: closure.add(u); changed=True
        m = mrev(w)
        if m not in closure: closure.add(m); changed=True
print("closure size:", len(closure), "(total words of len<=5: 1365)")
from collections import Counter
print(Counter(len(w) for w in closure))
