import subprocess
base='<project>'
out = subprocess.run(['grep','-n','mp.dps\\|ctx.prec\\|mp.prec\\|setprec','%s/paper/squircle_note/scripts/qsqrt2.py'%base,'%s/research/08_theory/A3_W4_proof/w4_right_coproduct.py'%base],capture_output=True,text=True)
print(out.stdout, out.stderr)
# head of w4_right_coproduct imports
out = subprocess.run(['sed','-n','1,40p','%s/research/08_theory/A3_W4_proof/w4_right_coproduct.py'%base],capture_output=True,text=True)
print(out.stdout)
