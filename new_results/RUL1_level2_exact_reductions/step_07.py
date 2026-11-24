import subprocess
print(subprocess.run(['bash','-c','cat "<project>/research/08_theory/K1f_hexagon_order5/result.md"'], capture_output=True, text=True).stdout)
