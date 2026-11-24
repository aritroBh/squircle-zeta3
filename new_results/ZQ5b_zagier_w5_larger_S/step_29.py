import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
dbg = r'''
from fractions import Fraction
import mpmath as mp
import sys
sys.path.insert(0,'.')
print("dps at import:", mp.mp.dps)
import w5blib as W
from w5blib import QD
eps = QD(Fraction(3,2),Fraction(1,2),13)
v = eps.emb(150)
print("direct emb(150) call at dps=15 default:")
mp.mp.dps=150
print("  vs 150-digit recompute:", mp.nstr(abs(v - abs(eps.emb(150))),5))
# now trace inside
old=mp.mp.dps; mp.mp.dps=150
w = mp.mpf(eps.a.numerator)/eps.a.denominator + (mp.mpf(eps.b.numerator)/eps.b.denominator)*mp.sqrt(eps.D)
print("  manual at 150:", mp.nstr(abs(w - (mp.mpf(3)+mp.sqrt(13))/2),5))
mp.mp.dps=old
'''
open('dbg5.py','w').write(dbg)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','dbg5.py'], capture_output=True, text=True, cwd='.')
print(r.stdout); print("ERR:", r.stderr[-800:])
# also grep the emb source as written
import subprocess
print(subprocess.run(['grep','-n','def emb','-A','5','w5blib.py'],capture_output=True,text=True,cwd='.').stdout)
