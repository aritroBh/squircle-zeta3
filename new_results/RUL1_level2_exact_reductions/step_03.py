import subprocess
print(subprocess.run(['bash','-c','cat "<project>/research/08_theory/D12_pi_p/result.md"'], capture_output=True, text=True).stdout)
