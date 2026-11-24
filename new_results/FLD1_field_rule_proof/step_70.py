import subprocess, time, json
WD='<project>/research/08_theory/FLD1_field_rule_proof'
VENV='<project>/.venv/bin/python'
# generate k=6 reference table (n<=60) by calling the mu() from step_21 without the N=100 verification block
code = r'''
import sys, json, time
src=open("step_21_fastmu.py").read()
ns={}
exec(src[:src.index("if __name__")], ns)
t=time.time()
KMAX=6; NMAX=60
D=ns['moments'](KMAX)
table={}
for k in range(KMAX+1):
    table[k]={}
    for n in range(1,NMAX+1):
        P=ns['mu'](k,n,KMAX,D)
        table[k][n]={','.join(map(str,e)):str(v) for e,v in P.items()}
print("computed k<=6 n<=60 in",round(time.time()-t,1),"s")
json.dump({str(k):{str(n):table[k][n] for n in table[k]} for k in table}, open('mu_exact_K6_N60.json','w'))
print("saved mu_exact_K6_N60.json")
print("mu_{6,2} =", table[6][2])
'''
open(WD+'/tmp_k6tab.py','w').write(code)
t=time.time()
r = subprocess.run(['nice',VENV,'tmp_k6tab.py'],capture_output=True,text=True,cwd=WD)
print("wall",round(time.time()-t,1),"s"); print(r.stdout); print(r.stderr[-800:] if r.returncode else '')