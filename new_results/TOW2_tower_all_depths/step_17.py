import shutil, os
base = '<project>/research/08_theory'
work = f'{base}/TOW2_tower_all_depths'
for f in ['towlib.py','fastred.py','genintegrand.py','buildintegrand.py','proof_common.py']:
    shutil.copy(f'{base}/TOW1_tower_law/{f}', f'{work}/{f}')
print(os.listdir(work))
# sanity: import and check data access
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, work)
os.chdir(work)
from towlib import *
print('keys degree 3:', keys_of_degree(3))
print('u(2,1, k=5, delta=2):', u(2,1,5,2))
print('rho((2,1),2):', rho((2,1),2), 'closed:', rho_c((2,1),2))
print('match:', rho((2,1),2)==rho_c((2,1),2))
open('progress.log','a').write('step_03: workspace set up, towlib copied, data access OK\n')
