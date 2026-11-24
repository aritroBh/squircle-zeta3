import subprocess
# Look at the paper section on corners of any angle
r = subprocess.run(['grep','-n','conj:field\|thm:flat\|thm:M2alpha\|prop:M3angle\|Corners of Any Angle\|Smoothed Corners','<project>/paper/ieee/main.tex'],capture_output=True,text=True)
print(r.stdout[:5000])