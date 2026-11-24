import sys, subprocess
dbg = r'''
import sys; sys.path.insert(0,'.')
import mpmath as mp
import w5lib as W
from w5lib import QD
C = W.letters_info(8)
eps = C['letters'][1]
x = eps**5
mp.mp.dps = 150
u = C['unit_emb']
print("unit_emb:", mp.nstr(u, 40))
v = y = x.emb(150)
print("x.emb:", mp.nstr(v, 40))
print("exact 41+29*sqrt2:", mp.nstr(41+29*mp.sqrt(2), 40))
l1 = mp.log(abs(v)); l2 = mp.log(u)
print("l1:", mp.nstr(l1, 40))
print("l2:", mp.nstr(l2, 40))
print("l2*5:", mp.nstr(l2*5, 40))
t = l1/l2
print("t:", mp.nstr(t, 40))
print("t-5:", mp.nstr(t-5, 10))
# check emb code path directly
old = mp.mp.dps; mp.mp.dps = 150
vv = mp.mpf(x.a.numerator)/x.a.denominator + (mp.mpf(x.b.numerator)/x.b.denominator)*mp.sqrt(x.D)
mp.mp.dps = old
print("manual emb:", mp.nstr(vv, 40))
print("diff emb-manual:", mp.nstr(v-vv, 10))
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',dbg],capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-1200:])
