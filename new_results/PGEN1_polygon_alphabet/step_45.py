# step_21: checks: (i) wedge3(delta)=0 for all polygons; (ii) 12-gon 7-point form?
import os, pickle
from fractions import Fraction as Fr
base = "<project>"
DOD = os.path.join(base,"research/08_theory/DOD5_dodecagon_termwise")
import subprocess
print(subprocess.run(["grep","-rl","seven\|7-pt\|7 pt\|7-point","--include=*.py",DOD],capture_output=True,text=True).stdout[:500])
print(subprocess.run(["grep","-rn","7 points\|seven points\|7-point","--include=*.md",DOD],capture_output=True,text=True).stdout[:800])
# check verify_dod5.py content
src = open(os.path.join(DOD,"verify_dod5.py")).read()
print(src[:3000])
