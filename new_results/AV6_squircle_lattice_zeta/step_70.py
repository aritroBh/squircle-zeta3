src = open("umach.py").read()
src = src.replace(
"""def eta(x):
    return (1-power(2, 1-x))*zeta(x)""",
"""def eta(x):
    from mpmath import log as mlog
    if abs(x-1) < mpf('1e-60'): return mlog(2)
    return (1-power(2, 1-x))*zeta(x)""")
open("umach.py","w").write(src)
test = r'''
from umach import *
from squircle import setdps
from mpmath import mp, nstr, zeta, gamma, power, mpf
import json
setdps(80)
s = 4
u0 = gamma(s-1)*zeta(s-1)*eta(s)
us = [u0] + [u_1_r(s, r) for r in range(1, 7)]
for r, u in enumerate(us):
    print(f"u_{{1,{r}}} = {nstr(u, 22)}")
# D(p) data
data = {int(k): v for k, v in json.load(open("Uq_data_s4.json")).items()}
print()
print("check D(p) = sum u_r p^{-r} vs data:")
for p_ in sorted(data.keys()):
    Ddata = mpf(data[p_][0])*power(p_, s-1)
    Dfit = sum(us[r]*power(p_, -r) for r in range(len(us)))
    print(f"p={p_}: data={nstr(Ddata, 18)}  series={nstr(Dfit, 18)}  diff={nstr(Ddata-Dfit, 5)}")
'''
exec(test)
