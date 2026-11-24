d='../AV5_layer_universality/'
import subprocess
out = subprocess.run(['grep','-n','pinf']+[d+f'step_{i}.py' for i in range(58,80)], capture_output=True, text=True)
print(out.stdout[:3000])
