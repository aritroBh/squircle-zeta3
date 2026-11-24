src = open("Uq_direct.py").read()
src = src.replace(
"""    I0q = fsum([((-1)**(n+1)) * _stirling_ratio(n, q) / (n*p+1) for n in range(q, 300)])""",
"""    I0q = fsum([_stirling_ratio(n, q) / (n*p+1) for n in range(q, 300)])""")
open("Uq_direct.py","w").write(src)
test = r'''
from squircle import setdps
from Uq_direct import Uq_direct
from mpmath import mp, nstr, mpf
setdps(60)
print("fixed U1(128):", nstr(Uq_direct(128, 4, 1, K=768, R=10), 25))
print("brute        : 0.000001035089172880316984234495")
print("fixed U1(20) :", nstr(Uq_direct(20, 4, 1, K=400), 25))
print("brute (old)  : 0.000212731477273906097966957952356")
'''
exec(test)
