import subprocess
print(subprocess.run(['bash','-c','sed -n 1,120p "<project>/research/08_theory/K1g_hexagon_order5_assembly/step_51.py"'], capture_output=True, text=True).stdout)
