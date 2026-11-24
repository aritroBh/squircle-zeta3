import os, shutil, time
base = '<project>/research/08_theory'
cwd = os.getcwd()
# find d5 candidates
for root, dirs, files in os.walk(f'{base}/ZQ5b_zagier_w5_larger_S'):
    if 'candidates.json' in files and 'd5' in root:
        print("FOUND d5:", root, files)
r = os.path.join(base,'ZQ5b_zagier_w5_larger_S','out')
print("out contents:", os.listdir(r))
