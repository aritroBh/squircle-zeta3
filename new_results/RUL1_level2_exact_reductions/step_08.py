import subprocess
print(subprocess.run(['bash','-c','cat "<project>/research/08_theory/K1g_hexagon_order5_assembly/result.md"'], capture_output=True, text=True).stdout)
