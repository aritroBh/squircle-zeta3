d='../AV5_layer_universality/'
import subprocess
out = subprocess.run(['grep','-n','hex3_epss\|oct4_epss\|def P_hex\|def P_oct\|q = \|q=\|1/(p-1)\|p/(p-1)']+[d+f'step_{i}.py' for i in range(86,118)], capture_output=True, text=True)
print(out.stdout[:4000])
