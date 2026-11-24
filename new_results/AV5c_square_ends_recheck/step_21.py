d='../AV5_layer_universality/'
import subprocess
# find the lame pinf steps and hexagon/octagon PSLQ steps
out = subprocess.run(['grep','-l','pslq'] + [d+f'step_{i}.py' for i in range(56,118)], capture_output=True, text=True)
print(out.stdout)
