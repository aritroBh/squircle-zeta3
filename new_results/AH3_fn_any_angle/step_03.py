
import os, re
os.chdir('<project>/research/08_theory/AH3_fn_any_angle')
with open('<project>/paper/ieee/main.tex') as fh:
    tex = fh.read()
# Save the two key sections to files for reading
open('sec_corners.tex','w').write(tex[83070:119070])
open('sec_fn.tex','w').write(tex[119070:121338])
open('sec_future.tex','w').write(tex[141910:147392])
print(tex[83070:88500])
