base = '<project>'
print('W4pp_1060 head:', open(base+'/research/08_theory/M4a_decomposition/W4pp_1060.txt').read()[:120])
print()
print('quad_90:', open(base+'/research/08_theory/M4a_decomposition/quad_90.json').read()[:200])
print()
import subprocess
# find Theorem W4 display in main.tex
out = subprocess.run(['grep','-n','Li_4','<project>/paper/squircle_note/main.tex'],capture_output=True,text=True).stdout
print(out[:2000])
