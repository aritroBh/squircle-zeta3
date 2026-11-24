import subprocess, sys
dbg = r'''
import sys; sys.path.insert(0,'.')
import mpmath as mp
import w5lib
from w5lib import QD
# monkeypatch a tracing emb
import mpmath as mp2
def emb_trace(s, prec=80):
    old = mp2.mp.dps; mp2.mp.dps = max(prec, 60)
    a = mp2.mpf(s.a.numerator)/s.a.denominator
    b = (mp2.mpf(s.b.numerator)/s.b.denominator)
    sq = mp2.sqrt(s.D)
    sys.stderr.write(f"dps={mp2.mp.dps} sq={mp2.nstr(sq,30)}\n")
    v = a + b*sq
    mp2.mp.dps = old
    return +v
x = QD(1,1,2)
mp.mp.dps = 15
v1 = emb_trace(x, 150)
mp.mp.dps = 50
print("traced emb:", mp.nstr(v1, 45))
v2 = x.emb(150)
print("method emb:", mp.nstr(v2, 45))
print("same object code?", QD.emb.__code__.co_filename, QD.emb.__code__.co_firstlineno)
print("mp in w5lib is mpmath?", w5lib.mp is mp, w5lib.mp.__name__, id(w5lib.mp), id(mp))
print("QD module file:", sys.modules['w5lib'].__file__)
import os
print("pycache:", os.listdir('__pycache__') if os.path.exists('__pycache__') else None)
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',dbg],capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-1200:])
