
import os
for f in ['D13b_Mk_constant.run','D13d_decay_theorem.run']:
    p = 'research/08_theory/'+f
    txt = open(p).read()
    print('='*20, f, len(txt))
    print(txt[:4500])
