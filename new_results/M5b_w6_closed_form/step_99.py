import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
src=open(os.path.join(m5b,'m5b_indep_check.py')).read()
src=src.replace("""def lchi5():
    return (mp.zeta(5, mp.mpf(1)/8) + mp.zeta(5, mp.mpf(3)/8)
            - mp.zeta(5, mp.mpf(5)/8) - mp.zeta(5, mp.mpf(7)/8))/8**5""",
"""def lchi5():
    # project convention: chi of Q(sqrt2), +1 on n = +/-1 mod 8, -1 on n = +/-3 mod 8
    return (mp.zeta(5, mp.mpf(1)/8) - mp.zeta(5, mp.mpf(3)/8)
            - mp.zeta(5, mp.mpf(5)/8) + mp.zeta(5, mp.mpf(7)/8))/8**5""")
open(os.path.join(m5b,'m5b_indep_check.py'),'w').write(src)
import subprocess, mpmath as mp, json
mp.mp.dps=60
L = (mp.zeta(5, mp.mpf(1)/8) - mp.zeta(5, mp.mpf(3)/8) - mp.zeta(5, mp.mpf(5)/8) + mp.zeta(5, mp.mpf(7)/8))/8**5
d=json.load(open(m5b+'/Lchi8_5_1800_planner.json'))
print('fixed lchi vs planner:', mp.nstr(abs(L-mp.mpf(d['value'][:62])),5))
print(open(m5b+'/w6_lll_1300.out').read()[-300:])
print('json1300:', os.path.exists(m5b+'/w6_lll_final_1300.json'))
