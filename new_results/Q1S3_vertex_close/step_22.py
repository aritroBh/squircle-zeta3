import sys
sys.path.insert(0, '.')
from q1s_exact import *
# self-tests
assert fnorm(p2_)==2
assert fnorm(q1_)==9 and fnorm(q2_)==9
assert fnorm(a_)==1
assert feq(q1_, fmul(i_, fsub(ONE, fmul(i_,r2))))
assert feq(q2_, fmul(fneg(i_), fadd(ONE, fmul(i_,r2))))
for nm,x in [('a',a_),('p2',p2_),('q1',q1_),('q2',q2_),('sqrt2',r2),('i',i_),('1+i',fadd(ONE,i_))]:
    print(nm, factor_sunit(x))
print("exact engine OK")
