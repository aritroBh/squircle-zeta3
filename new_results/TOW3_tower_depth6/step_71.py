import sys
base = '<project>/research/08_theory'
wts = f'{base}/WTS_weight_k_structure/theory'
sys.path.insert(0, wts)
from fastred import DL1, DL2
def show(P):
    print('real:', P[0]); print('ipos:', P[1]); print('ineg:', P[2]); print('const:', P[3])
print('DL1:'); show(DL1)
print('DL2:'); show(DL2)
import subprocess
print(subprocess.run(['grep','-n','DL1\s*=\|DL2\s*=', f'{wts}/fastred.py'],capture_output=True,text=True).stdout)
src = open(f'{wts}/fastred.py').read()
i = src.find('DL1')
print(src[i-200:i+400])
