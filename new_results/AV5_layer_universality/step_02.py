import os, shutil
os.chdir('<project>/research/08_theory/AV5_layer_universality')
# remove the accidental nested dir
if os.path.isdir('research'):
    shutil.rmtree('research')
print(os.listdir('.'))
if os.path.exists('step_01.py'):
    print(open('step_01.py').read()[:3000])
