import subprocess
out = subprocess.run(['grep','-l','sol_B0','../OCT9b_octagon_debug/step_*.py'],capture_output=True,text=True)
print(out.stdout)
out = subprocess.run(['grep','-l','sol_B123','../OCT9b_octagon_debug/step_*.py'],capture_output=True,text=True)
print(out.stdout)
out = subprocess.run(['grep','-ln','W_parts_e0','../OCT9b_octagon_debug/step_*.py'],capture_output=True,text=True)
print(out.stdout)
