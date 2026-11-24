import subprocess
print(subprocess.run(['bash','-c','sed -n 1,90p "<project>/research/08_theory/K1g_hexagon_order5_assembly/step_50.py"; echo ====ATOM_EVAL; sed -n 1,80p "<project>/research/08_theory/K1g_hexagon_order5_assembly/atom_eval.py"'], capture_output=True, text=True).stdout)
