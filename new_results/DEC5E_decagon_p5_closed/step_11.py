import subprocess
out = subprocess.run(['grep','-rln','wparts','../OCT5D3_octagon_p5_closed/'],capture_output=True,text=True)
print(out.stdout)
out = subprocess.run(['grep','-rln','Wt','../OCT5D3_octagon_p5_closed/step_*.py'],capture_output=True,text=True)
print(out.stdout[:2000])
