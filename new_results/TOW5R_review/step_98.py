import pickle, sympy as sp, sys
wd = '<project>/research/08_theory/TOW5R_review'
d, m = sp.symbols('d m')
raw = pickle.load(open(wd+'/tow5_copy/forced_F910.pkl','rb'))
for key in [(9,1),(10,1),(10,2)]:
    v = sp.sympify(raw[key], locals={'d':d,'m':m,'I':sp.I})
    n, dn = sp.fraction(sp.together(v))
    nc, dc = sp.fraction(sp.cancel(v))
    print(key, 'stored-str terms in Add:', len(sp.Add.make_args(sp.expand(n))) if True else '?',
          '| together-num terms:', len(sp.Poly(n, d, m).terms()),
          '| cancel-num terms:', len(sp.Poly(nc, d, m).terms()),
          '| den:', sp.factor(dn))
