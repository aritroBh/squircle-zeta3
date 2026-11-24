
import os, subprocess
os.chdir('<project>/research/08_theory/AH3_fn_any_angle')
files = sorted([f for f in os.listdir('.') if f.startswith('step_') and f.endswith('.py')], key=lambda s: int(s[5:7]))
print(files)
print(subprocess.run(['wc','-l']+files, capture_output=True, text=True).stdout)
