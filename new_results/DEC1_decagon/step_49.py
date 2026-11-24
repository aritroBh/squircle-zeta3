import subprocess
d='research/08_theory/OCT6_octagon_closed/'
# find the step with the analytic proof identities (Landen-sum, inversion, difference)
out=subprocess.run(['grep','-ln','DIFF\|Landen\|difference identity',*[d+f'step_{i:02d}.py' for i in range(6,64)]],capture_output=True,text=True)
print(out.stdout)
out=subprocess.run(['grep','-n','def R_L2\|def R_L1\|def DIFF\|R_L2 =\|R_L1 =\|DIFF =\|Landen-sum\|inversion',*[d+f'step_{i:02d}.py' for i in range(6,64)]],capture_output=True,text=True)
print(out.stdout[:4000])