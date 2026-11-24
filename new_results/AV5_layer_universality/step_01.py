import os, subprocess
print(os.getcwd())
print(os.listdir('.'))
os.makedirs('research/08_theory/AV5_layer_universality', exist_ok=True)
os.chdir('research/08_theory/AV5_layer_universality')
print(os.getcwd())
# check what's around
print(subprocess.run(['ls','../../../..'],capture_output=True,text=True).stdout[:2000])
with open('progress.log','a') as f:
    f.write('step_00: setup, cwd initialized\n')
