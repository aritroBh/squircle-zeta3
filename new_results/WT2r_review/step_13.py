import subprocess
base = '<project>/research/08_theory/WT2b_stabilization'
r = subprocess.run(['grep', '-l', 'record_W2b', f'{base}/step_0*.py', f'{base}/step_1*.py', f'{base}/step_2*.py'], capture_output=True, text=True)
print(r.stdout)
r2 = subprocess.run(['grep', '-rn', "record_W2b.pkl','wb'", base], capture_output=True, text=True)
print("WRITER:", r2.stdout)
r3 = subprocess.run(['grep', '-rn', "record_W2b.pkl\",'wb'", base], capture_output=True, text=True)
print("WRITER2:", r3.stdout)
