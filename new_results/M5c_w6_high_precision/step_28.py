import json
sol=json.load(open('w6_symbol_solve.json'))
print('atoms:',len(sol['names']),'rank:',sol['rank'],'kernel dim:',len(sol['Ksat']))
print('li6_args:',sol['li6_args'])
import itertools
from fractions import Fraction as F
# light check: xp/xq are rationals, kernel rows are ints of len 149
assert all(len(r)==149 for r in sol['Ksat'])
assert all(isinstance(v,int) for r in sol['Ksat'] for v in r)
print('xp nonzero:',sum(1 for c in sol['xp'] if F(c)!=0),'xq nonzero:',sum(1 for c in sol['xq'] if F(c)!=0))
mx=max(max(abs(F(c).numerator),1) for c in sol['xp']+sol['xq'])
print('particular max|num|:',mx)
