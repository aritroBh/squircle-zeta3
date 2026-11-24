import os, json
os.chdir('<project>')
print(os.getcwd())
base='research/08_theory/M4e_w5_lchi'
for root,dirs,files in os.walk(base):
    for f in files:
        p=os.path.join(root,f)
        print(p, os.path.getsize(p))
