import subprocess
t5 = '<project>/research/08_theory/TOW5_tower_depth10'
# 67/67 claim source
s68 = open(t5+'/step_68.py').read()
print(s68[:600])
print('...')
r = subprocess.run(['grep','-c','67/67', t5+'/transcript.json'], capture_output=True, text=True)
print('67/67 in transcript:', r.stdout)
r = subprocess.run(['grep','-n','fitted+verified', t5+'/tow2lib.py'], capture_output=True, text=True)
print('stale docstring:', r.stdout)
