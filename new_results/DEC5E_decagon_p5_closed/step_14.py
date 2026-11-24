import subprocess
out = subprocess.run(['grep','-ln','wparts_true','../OCT5D3_octagon_p5_closed/step_*.py'],capture_output=True,text=True)
print(out.stdout)
print(open('../OCT5D3_octagon_p5_closed/step_53.py').read()[:4000])