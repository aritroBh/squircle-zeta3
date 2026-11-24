import sys, subprocess
dbg = r'''
import sys; sys.path.insert(0,'.')
import mpmath as mp
from w5lib import QD
import w5lib as W
x = QD(1,1,2)
v = x.emb(150)
mp.mp.dps = 60
print("emb(150):", mp.nstr(v, 55))
mp.mp.dps = 150
print("1+sqrt2 :", mp.nstr(1+mp.sqrt(2), 55))
mp.mp.dps = 60
# now what letters_info returns
C = W.letters_info(8)
print("unit   :", C['unit'], " D=", C['unit'].D)
print("unit_emb:", mp.nstr(C['unit_emb'], 55))
# reproduce emb manually at dps=15 outer
mp.mp.dps = 15
v2 = x.emb(150)
mp.mp.dps = 60
print("emb(150) with outer dps 15:", mp.nstr(v2, 55))
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',dbg],capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-1200:])
