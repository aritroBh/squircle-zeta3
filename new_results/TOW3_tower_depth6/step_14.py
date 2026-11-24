base = '<project>/research/08_theory'
t2 = f'{base}/TOW2_tower_all_depths'
import subprocess
out = subprocess.run(['grep','-rl','towers_D8\|TOWN1\|\.pkl', t2], capture_output=True, text=True)
print(out.stdout[:3000])
# how did TOW2 load exact tower data?
out2 = subprocess.run(['grep','-rn','H\[' , t2+'/step_77.py'], capture_output=True, text=True)
print(out2.stdout[:2000])
