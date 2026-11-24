import sys, json
sys.path.insert(0, ".")
from fractions import Fraction as Fr
from cor5_sym import W4Quotient, li_sym

Q = W4Quotient(2)
sols = json.load(open("sunit23_solutions.json"))
V23 = {Fr(xs): (Fr(fx[1]), Fr(fx[2])) for xs, fx, fy in sols}
OM23 = {Fr(xs): (Fr(fy[1]), Fr(fy[2])) for xs, fx, fy in sols}
def li4v(x): return Q.coords(li_sym(4, V23[x], OM23[x]))

# numeric F4 relation's Li4 part
num_li4 = {Fr(1,3): Fr(103,8), Fr(2,3): Fr(-29,2), Fr(3,4): Fr(-21,4),
           Fr(1,2): Fr(-31,6), Fr(1,9): Fr(-45,64)}
tot = [Fr(0)]*3
for x, c in num_li4.items():
    v = li4v(x)
    for i in range(3):
        tot[i] += c * v[i]
print("projection of numeric Li4 part:", tuple(tot))
print("vF (from symbol):                 ", (Fr(5,2), Fr(1,2), Fr(0)))
# the symbol-proved combination
sym_li4 = {Fr(1,3): Fr(-1,2), Fr(2,3): Fr(-3,2), Fr(3,4): Fr(-1,2)}
tot2 = [Fr(0)]*3
for x, c in sym_li4.items():
    v = li4v(x)
    for i in range(3):
        tot2[i] += c * v[i]
print("projection of symbol combination: ", tuple(tot2))
# difference: must be a Li4 functional-equation (zero projection) combination
diff = {x: num_li4.get(x, Fr(0)) - sym_li4.get(x, Fr(0)) for x in set(num_li4) | set(sym_li4)}
totd = [Fr(0)]*3
for x, c in diff.items():
    if c:
        v = li4v(x)
        for i in range(3):
            totd[i] += c * v[i]
print("projection of (numeric - symbol): ", tuple(totd))
