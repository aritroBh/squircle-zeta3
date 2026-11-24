
wd = 'research/08_theory/K2b_corner_phase/'
s = open(wd+'k2utils.py').read()
s = s.replace("import re, mpmath\n", "import re\nimport mpmath\nfrom mpmath import mp\n")
open(wd+'k2utils.py','w').write(s)
import subprocess
print(open(wd+'k2utils.py').read()[:200])
