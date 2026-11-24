import subprocess
t5 = '<project>/research/08_theory/TOW5_tower_depth10'
base = '<project>/research/08_theory'
r = subprocess.run(['cmp', '-s', t5+'/trunc_engine.py', base+'/TOWN1_tower_large_n/trunc_engine.py'])
print('trunc_engine SAME as TOWN1:', r.returncode == 0)
r = subprocess.run(['cmp', '-s', t5+'/driver.py', base+'/TOWN1_tower_large_n/driver.py'])
print('driver SAME:', r.returncode == 0)
# find check_data.py source
import os
for d in ['TOW3R_review','TOW3_tower_depth6','TOW4_tower_depth8','TOW4R_review']:
    p = base+'/'+d+'/check_data.py'
    if os.path.exists(p):
        r = subprocess.run(['cmp','-s',t5+'/check_data.py',p])
        print('check_data SAME as', d, ':', r.returncode == 0)
