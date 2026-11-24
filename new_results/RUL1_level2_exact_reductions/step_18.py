import subprocess
print(subprocess.run(['bash','-c','ls "<project>/research/08_theory/D12_pi_p/"; sed -n 1,150p "<project>/research/08_theory/D12_pi_p/step_04.py"'], capture_output=True, text=True).stdout)
