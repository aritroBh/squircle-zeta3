import os, pickle, sympy as sp
os.chdir('research/08_theory/DEC1_decagon')
exact = pickle.load(open('E3n_exact_pq.pkl','rb'))
r5 = sp.sqrt(5)
tt, L0,L1,L2,L3,L4 = sp.symbols('t L0 L1 L2 L3 L4')
lt,l1t,w1s,w4s,l2,lf = sp.symbols('lt l1t w1 w4 l2 lf')
Ls=[L0,L1,L2,L3,L4]
E = sp.Integer(0)
for key,row in exact.items():
    mon = sp.Integer(1)
    for i,e in enumerate(key): mon *= Ls[i]**e
    for nn in range(4):
        p,q = row[nn]
        if p==0 and q==0: continue
        E += (p+q*r5)*tt**nn*mon
E = sp.expand(E)
print("terms in L-basis:", len(sp.Add.make_args(E)))
sub = {L0: lt-l2-2*lf, L1: l1t-l2-2*lf, L2: w1s-l2-2*lf, L3: w4s-l2-lf, L4: -l2-lf}
E2 = sp.expand(E.xreplace(sub))
print("terms in reduced basis:", len(sp.Add.make_args(E2)))
mons = {}
for term in sp.Add.make_args(E2):
    # decompose: coeff(p+qr5) * tt^n * integ-monomial(lt,l1t,w1s,w4s) * letters(l2,lf)
    d = term.as_powers_dict()
    n_ = int(d.get(tt,0))
    abc = (int(d.get(lt,0)), int(d.get(l1t,0)), int(d.get(w1s,0)), int(d.get(w4s,0)))
    letters = (int(d.get(l2,0)), int(d.get(lf,0)))
    rest = sp.simplify(term/(tt**n_ * lt**abc[0]*l1t**abc[1]*w1s**abc[2]*w4s**abc[3] * l2**letters[0]*lf**letters[1]))
    key = (n_, abc, letters)
    mons[key] = sp.expand(mons.get(key,0)+rest)
print("distinct (n, integ, letter) monomials:", len(mons))
# classify
from collections import Counter
cls = Counter()
for (n_,abc,letters) in mons:
    a,b,c1,c4 = abc
    if c1==0 and c4==0: cls['local']+=1
    elif c1>0 and c4==0: cls['single_x1']+=1
    elif c4>0 and c1==0: cls['single_x4']+=1
    else: cls['mixed']+=1
print(cls)
# mixed types
mt = sorted(set((abc[0],abc[1],abc[2],abc[3]) for (n_,abc,letters) in mons if abc[2]>0 and abc[3]>0))
print("mixed integral types (a,b,c1,c4):", mt)
ns_needed = {}
for mtype in mt:
    ns_needed[mtype] = sorted(set(n_ for (n_,abc,letters) in mons if abc==mtype))
print("n values per mixed type:", ns_needed)
pickle.dump(mons, open('mons.pkl','wb'))
with open('progress.log','a') as f:
    f.write(f"step: basis change done; {len(mons)} monomials; classes {dict(cls)}; mixed types {mt}; n per type {ns_needed}\n")