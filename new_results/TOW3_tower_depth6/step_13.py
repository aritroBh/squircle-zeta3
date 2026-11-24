base = '<project>/research/08_theory'
t2 = f'{base}/TOW2_tower_all_depths'
import subprocess
# find which steps contain rho verification vs data, and tower law check
out = subprocess.run(['grep','-l','towers_D8_K52'] + [f'{t2}/step_{n:02d}.py' for n in range(1,87)], capture_output=True, text=True)
print(out.stdout)
