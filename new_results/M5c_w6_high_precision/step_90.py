import sys, pickle, json
sys.dont_write_bytecode=True
work='<project>/research/08_theory/M5c_w6_high_precision'
sys.path.insert(0,work)
import mpmath as mp
mp.mp.dps=3620
from m5b_common import pool_zs
from w6_solve import atoms_w6, BASE_LI6
pool=pool_zs()
li6=list(BASE_LI6)+[pool[i] for i in [14,15,13,9,7,6]]
basis=atoms_w6(li6)
n2v={b[0]: mp.mpf(v) for b,v in zip(basis, pickle.loads(open(work+'/w6_atomvals_3500_x14_15_13_9_7_6.pkl','rb').read()))}
def L(z): return n2v[f'Li6({z})']
import math
# Li6(b^2) vs 32(Li6(b)+Li6(-b))
d1 = L('(3-2*sqrt2)') - 32*(L('(-1+1*sqrt2)')+L('(1-1*sqrt2)'))
# Li6(b^4) vs 32(Li6(b^2)+Li6(-b^2)) and vs 1024 base combo
d2 = L('(17-12*sqrt2)') - 32*(L('(3-2*sqrt2)')+L('(-3+2*sqrt2)'))
d2b= L('(17-12*sqrt2)') - (1024*L('(-1+1*sqrt2)')+1024*L('(1-1*sqrt2)')+32*L('(-3+2*sqrt2)'))
# Li6(-1/sqrt2) vs Li6(1/2)/32 - Li6(1/sqrt2)
d3 = L('-1/2*sqrt2') - (L('1/2')/32 - L('1/2*sqrt2'))
print('Li6(b^2)  = 32(Li6(b)+Li6(-b))        residual:', mp.nstr(abs(d1),5))
print('Li6(b^4)  = 32(Li6(b^2)+Li6(-b^2))    residual:', mp.nstr(abs(d2),5))
print('Li6(b^4)  = 1024(base combo)          residual:', mp.nstr(abs(d2b),5))
print('Li6(-1/r2)= Li6(1/2)/32 - Li6(1/r2)   residual:', mp.nstr(abs(d3),5))
