import subprocess
print(subprocess.run(['bash','-c','cat "<project>/research/08_theory/E5e_pslq_340/result.md"'], capture_output=True, text=True).stdout)
