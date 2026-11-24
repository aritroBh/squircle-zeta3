src = open('step_10.py').read()
old = "d5 = omreg_up**2*up(invp2reg)/2"
new = "d5 = omreg_up**2*up(invp2reg)/2*(y_J*omreg_up**2*up(invp2reg)/2).exp()  # |e^x-1|<=|x|e^|x|, x=om^2 invp2/2"
assert old in src; src = src.replace(old, new)
open('step_10.py','w').write(src)
# regenerate J192 variant
v = src.replace('J = 48; zstar = arb("3.0")', 'J = 192; zstar = arb("3.0")')
v = v.replace('open("final_cell.json","w")', 'open("final_cell_J192.json","w")')
open('step_10_J192.py','w').write(v)
import subprocess, sys
for f in ['step_10.py', 'step_10_J192.py']:
    r = subprocess.run([sys.executable, f], capture_output=True, text=True, timeout=600)
    for line in r.stdout.splitlines():
        if 'FINAL CELL' in line or 'd1..d6' in line or 'Dreg' in line:
            print(f"{f}:", line)
    if r.stderr: print(f, "STDERR:", r.stderr[-500:])
with open('progress.log','a') as f2:
    f2.write("turn9: d5 fix (|e^x-1|<=|x|e^|x|) applied; final cell rerun.\n")
