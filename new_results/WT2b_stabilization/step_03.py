import os
base = '<project>/research/08_theory'
print(os.getcwd())
print(os.listdir(base+'/WT2b_stabilization') if os.path.exists(base+'/WT2b_stabilization') else 'MISSING')
