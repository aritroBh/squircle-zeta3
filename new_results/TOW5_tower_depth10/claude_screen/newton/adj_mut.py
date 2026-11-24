# Mutation + independent integer-p evaluation for prove_adj910_fast.py
import sys, time
sys.dont_write_bytecode=True
W="<project>/research/08_theory/TOW5_tower_depth10"
sys.path.insert(0,W)
import sympy as sp
import tow5lib as T5
import prove_adj910_fast as P
p=P.p
def pointcheck(delta,s,pts):
    lhs=P.eR(delta,s,P.d,p)+P.eR(delta,s+1,P.d,p)
    terms=P.feed_terms(delta,s)
    out=[]
    for pv in pts:
        sub={p:sp.Integer(pv)}
        val=sp.nsimplify(sp.expand(sp.simplify(lhs.subs(sub)-sum(t.subs(sub) for t in terms))))
        out.append((pv,val))
    return out
mode=sys.argv[1]; delta=int(sys.argv[2]); s=int(sys.argv[3])
if mode=="base":
    print("fast method:"); 
    try: P.main(delta,s)
    except AssertionError: print("FAST METHOD ASSERTION FAILED (unexpected for base)")
    print("integer p direct substitution (lhs - sum terms):",pointcheck(delta,s,[4,5,6,7,9,12]))
elif mode=="tau":   # perturb tau2(2) by 1/256
    orig=T5.tau2f
    T5.tau2f=lambda n: orig(n)+(sp.Rational(1,256) if n==2 else 0)
    try: P.main(delta,s); print("MUTANT PASSED (vacuous!)")
    except AssertionError: print("mutant tau2(2)+1/256: fast method FAILS as required")
elif mode=="ehat": # perturb the cached ehat fraction
    orig=P.ehat
    P._ehat.clear()
    def bad(delta_,s_,conj=False):
        v=orig(delta_,s_,conj)
        return sp.cancel(v+(sp.Rational(1,7)/(P.p+1) if (delta_,s_,conj)==(4,1,False) else 0))
    P.ehat=bad; P._ed.clear()
    try: P.main(delta,s); print("MUTANT PASSED (vacuous!)")
    except AssertionError: print("mutant ehat(4,1)+1/(7(p+1)): fast method FAILS as required")
