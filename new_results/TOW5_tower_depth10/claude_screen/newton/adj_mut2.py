# targeted mutants for prove_adj910_fast.main(delta,s): perturb ONE input that the identity actually uses
import sys
sys.dont_write_bytecode=True
W="<project>/research/08_theory/TOW5_tower_depth10"
sys.path.insert(0,W)
import sympy as sp
import tow5lib as T5
import prove_adj910_fast as P
mode=sys.argv[1]; delta=int(sys.argv[2]); s=int(sys.argv[3]); a=int(sys.argv[4]); b=int(sys.argv[5]) if len(sys.argv)>5 else 0
eps=sp.Rational(1,256)
if mode=="tau1":
    o=T5.tau1f; T5.tau1f=lambda n: o(n)+(eps if n==a else 0)
elif mode=="tau2":
    o=T5.tau2f; T5.tau2f=lambda n: o(n)+(eps if n==a else 0)
elif mode=="ehat":   # ehat(a,b,conj=False) += 1/(7(p+1))
    o=P.ehat
    def bad(dl,ss,conj=False):
        v=o(dl,ss,conj)
        return sp.cancel(v+(sp.Rational(1,7)/(P.p+1) if (dl,ss,conj)==(a,b,False) else 0))
    P.ehat=bad; P._ed.clear()
try:
    P.main(delta,s); print(f"MUTANT {mode} {a} {b} at ({delta},{s}): PASSED (insensitive)")
except AssertionError:
    print(f"MUTANT {mode} {a} {b} at ({delta},{s}): FAILS as required")
