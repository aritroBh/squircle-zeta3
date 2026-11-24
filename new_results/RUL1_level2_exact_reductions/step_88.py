import subprocess
PY = '<project>/.venv/bin/python'
print(subprocess.run(['bash','-c','sed -n 1,90p "<project>/research/08_theory/K1g_hexagon_order5_assembly/step_47.py"'], capture_output=True, text=True).stdout)
