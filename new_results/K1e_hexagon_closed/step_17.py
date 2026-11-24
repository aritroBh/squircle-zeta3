
import subprocess, os
base='<project>'
out=subprocess.run(['grep','-rl','t3_hexagon','--include=*.py','.'],capture_output=True,text=True,cwd=base)
print(out.stdout[:2000], out.stderr[:500])
out=subprocess.run(['grep','-rl','P3_arc','--include=*.py','.'],capture_output=True,text=True,cwd=base)
print(out.stdout[:2000])
