import pickle, sympy as sp
wd = '<project>/research/08_theory/K1f_hexagon_order5'
jtab = pickle.load(open(f'{wd}/jtab5.pkl','rb'))
ATOMSYMS = sp.symbols('H1 Li2 Li3 H011 Li4 H0011 H0111 Li5 H00011 H00101 H00111 H01011 H01111')
z2,z3,z5 = sp.symbols('z2 z3 z5')
names = {str(s): s for s in ATOMSYMS}
loc = {str(s): s for s in ATOMSYMS}; loc.update({'z2':z2,'z3':z3,'z5':z5}); loc['x']=sp.Symbol('x')
used = {}
maxpole = 0
for k, v in jtab.items():
    e = sp.sympify(v, locals=loc)
    for term in sp.Add.make_args(e):
        pd = term.as_powers_dict()
        xx = pd.get(sp.Symbol('x'), 0)
        maxpole = min(maxpole, xx)
        for s in ATOMSYMS + (z2,z3,z5):
            if pd.get(s,0):
                used.setdefault(str(s), set()).add(k)
print("atoms appearing in the 210-key J-table (with # of keys):")
for s in ['H1','Li2','Li3','H011','Li4','H0011','H0111','Li5','H00011','H00101','H00111','H01011','H01111','z2','z3','z5']:
    if s in used:
        ks = used[s]
        print(f"  {s:8s}: {len(ks):3d} keys")
    else:
        print(f"  {s:8s}: ABSENT")
print("most negative x-power in J-table:", maxpole)
# also: do z5 or products appear? check for products of atoms (weight<=5 monomials)
import collections
monset = collections.Counter()
for k, v in jtab.items():
    e = sp.sympify(v, locals=loc)
    for term in sp.Add.make_args(e):
        pd = term.as_powers_dict()
        key = tuple(sorted((str(s), pd.get(s,0)) for s in ATOMSYMS+(z2,z3,z5) if pd.get(s,0)))
        monset[key] += 1
print("distinct atom monomials in J-table:", len(monset))
