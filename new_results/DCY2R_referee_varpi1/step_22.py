import subprocess
# CX6 proof: find the ray representation and rest bound
r = subprocess.run(['grep','-n','main\|10,600\|10600000\|rest\|Im.*tilde\|tilde I\|2/\\\\pi\|eq. (20)\|aR^2','../CX6_decay_rest/proof.md'],capture_output=True,text=True)
print(r.stdout[:6000])
print("="*60)
print(subprocess.run(['head','-5','../D13f_sweep/mk_deg9.tsv'],capture_output=True,text=True).stdout)
print(subprocess.run(['bash','-c','sed -n "155,160p" ../D13f_sweep/mk_deg9.tsv'],capture_output=True,text=True).stdout)
