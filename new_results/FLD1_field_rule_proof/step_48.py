import subprocess, os
# Look for G2_M4alpha
for root in ['<project>/research/10_frontier/','<project>/research/']:
    r = subprocess.run(['bash','-c',f'ls -d {root}*G2* {root}*M4* 2>/dev/null; find {root} -maxdepth 2 -iname "*g2*m4*" -o -iname "*m4alpha*" 2>/dev/null | head'],capture_output=True,text=True)
    if r.stdout.strip(): print(root, "->", r.stdout)
# also grep the paper for S-unit statement
r = subprocess.run(['grep','-n','S-unit\|S-units\|{2,3}-unit\|unit equation','<project>/paper/ieee/main.tex'],capture_output=True,text=True)
print("paper S-unit mentions:", r.stdout[:2000])