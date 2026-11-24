import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
dbg = r'''
import sys
from fractions import Fraction
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
F = W.make_field(13, split=[3], inert=[2])
ue = F['unit_emb']
print("unit_emb type:", type(ue))
mp.mp.dps = 150
direct = abs(F['unit'].emb(150))
print("diff unit_emb vs fresh:", mp.nstr(abs(ue - direct), 5))
print("prec of unit_emb:", mp.prec(ue) if hasattr(ue,'_prec') else '?')
y = QD(Fraction(119,2),Fraction(33,2),13)
ye = abs(y.emb(150))
print("log(y.emb):", mp.nstr(mp.log(ye), 30))
print("log(unit_emb):", mp.nstr(mp.log(ue), 30))
print("ratio:", mp.nstr(mp.log(ye)/mp.log(ue), 30))
print("log(direct):", mp.nstr(mp.log(direct), 30))
'''
open('dbg4.py','w').write(dbg)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','dbg4.py'], capture_output=True, text=True, cwd='.')
print(r.stdout); print("ERR:", r.stderr[-1500:])
