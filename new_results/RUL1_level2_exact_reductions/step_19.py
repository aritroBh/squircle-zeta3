import subprocess
print(subprocess.run(['bash','-c','sed -n 1,80p "<project>/research/08_theory/D12_pi_p/step_05.py"'], capture_output=True, text=True).stdout)
