import sys, subprocess
dbg = r'''
import sys; sys.path.insert(0,'.')
import mpmath as mp
from w5lib import QD
x = QD(1,1,2)
print("fields:", repr(x.a), repr(x.b), repr(x.D), type(x.a.numerator))
mp.mp.dps = 15
old = mp.mp.dps
mp.mp.dps = max(150,60)
print("dps now:", mp.mp.dps)
t1 = mp.mpf(x.a.numerator)/x.a.denominator
t2 = (mp.mpf(x.b.numerator)/x.b.denominator)*mp.sqrt(x.D)
print("t1:", mp.nstr(t1, 50))
print("sqrt(2):", mp.nstr(mp.sqrt(2), 50))
print("t2:", mp.nstr(t2, 50))
v = t1 + t2
print("v:", mp.nstr(v, 50))
mp.mp.dps = old
# compare with method result
w = x.emb(150)
mp.mp.dps = 60
print("method v:", mp.nstr(w, 50))
import inspect
import w5lib
print(inspect.getsource(w5lib.QD.emb))
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',dbg],capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-1200:])
