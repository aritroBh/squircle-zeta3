import subprocess
out = subprocess.run(['grep','-ln','wparts','../OCT5D3_octagon_p5_closed/step_*.py'],capture_output=True,text=True)
print(out.stdout)
