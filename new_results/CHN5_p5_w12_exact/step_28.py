import subprocess
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
# structure of R_factors.txt: headers only (degrees)
out=subprocess.run(['grep','-n','factor\|deg'], capture_output=True, text=True, input=open(W+'/R_factors.txt').read()).stdout
print(out[:3000])
print("total lines:", len(open(W+'/R_factors.txt').readlines()), "size:", len(open(W+'/R_factors.txt').read()))
