import sys, subprocess
dbg = r'''
import sys; sys.path.insert(0,'.')
import mpmath as mp
import w5lib as W
from w5lib import QD
C = W.letters_info(8)
w = C['letters'][0]; eps = C['letters'][1]
print("letters:", w, eps, "unit_idx:", C['unit_idx'])
x = eps**5
print("x =", x, "N =", x.norm())
print("pexps:", C['pexps'](x.norm()))
y = x / (w**0)
print("y =", y, "norm:", y.norm())
mp.mp.dps = 150
t = mp.log(abs(y.emb(150)))/mp.log(C['unit_emb'])
print("t =", mp.nstr(t, 30), "nint:", int(mp.nint(t)))
k = int(mp.nint(t))
print("y == eps**k:", y == eps**k, " y == -eps**k:", y == -(eps**k))
print("eps**k =", eps**k)
print("factor_sunit(x):", W.factor_sunit(x, C))
print("factor_sunit(-x):", W.factor_sunit(-x, C))
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',dbg],capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
