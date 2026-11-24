import sys, pickle, json
sys.dont_write_bytecode=True
work='<project>/research/08_theory/M5c_w6_high_precision'
sys.path.insert(0,work)
import mpmath as mp
mp.mp.dps=3620
from m5b_common import qnum
from w6_solve import n6_terms
from w6_values import leaf_recipes
recs=leaf_recipes()
saved=pickle.loads(open(work+'/w6_leaves_3500.pkl','rb').read())
leaves={r: mp.mpf(v) for r,v in zip(recs,saved)}
n6v=mp.mpf(0)
for m,n,x,y,c in n6_terms():
    n6v += qnum(c)*leaves[('mpl',m,n,x,y)]
n6ref=json.load(open('<project>/research/08_theory/M5a_next_coefficient/n6_values.json'))
w6=mp.mpf(json.load(open(work+'/w6_11800.json'))['value_full'])
mp.mp.dps=80
print('my N6   :', mp.nstr(n6v,40))
print('M5a N6  :', n6ref['N6'][:44])
print('diff    :', mp.nstr(abs(n6v-mp.mpf(n6ref['N6'])),5))
print('my W6-N6:', mp.nstr(w6-n6v,40))
print('M5a W-N6:', n6ref['W6_minus_N6'][:44])
print('diff    :', mp.nstr(abs((w6-n6v)-mp.mpf(n6ref['W6_minus_N6'])),5))
