import sys, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
Fraction = Fr
P = G.golden_point('gol1')
F = P.F
def nfadd(a, b): return (a[0]+b[0], a[1]+b[1])
def nfmul(a, b): return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
Z2 = (Fr(0), Fr(0))
RHO3 = (Fr(-2), Fr(1))
P5 = [((Fr(-233, 5), Fr(54, 5)), (Fr(3, 4), Fr(-1, 4))),
      ((Fr(-23, 5), Fr(-6, 5)), (Fr(1, 4), Fr(1, 4))),
      ((Fr(-107, 45), Fr(233, 135)), RHO3),
      ((Fr(-32, 5), Fr(-4, 5)), (Fr(3, 1), Fr(-1, 1)))]
SP5 = {}
for c, z in P5:
    for w, cf in G.li5_sym_nf(P, z).items():
        SP5[w] = nfadd(SP5.get(w, Z2), nfmul(c, cf))
cl_P5 = G.quotient5(G.ft_expand_sym_nf(SP5, 2))
print("class(P5) =", [G.nf_str(q) for q in cl_P5])
# target from scan: 
T = ((Fr(192), Fr(-24)), (Fr(0), Fr(0)), (Fr(-108), Fr(0)), (Fr(0), Fr(0)),
     (Fr(-147), Fr(-3)), (Fr(636), Fr(-69)))
print("target    =", [G.nf_str(q) for q in T])
# check against li5_classes_golden.json classes directly
li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
tot = [Z2]*6
for c, z in P5:
    clz = li5cl[z]
    tot = [nfadd(t, nfmul(c, q)) for t, q in zip(tot, clz)]
print("class via stored li5cl =", [G.nf_str(q) for q in tot])
