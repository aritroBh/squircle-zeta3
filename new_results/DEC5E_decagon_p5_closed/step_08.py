import subprocess
# Which DEC5D step computed P5/L from rewrite? grep
out = subprocess.run(['grep','-l','rewrite','../DEC5D_decagon_order5_depth/step_5*.py'],capture_output=True,text=True)
print(out.stdout)
out = subprocess.run(['grep','-rn','P5','../DEC5D_decagon_order5_depth/step_53.py'],capture_output=True,text=True)
print(out.stdout)
print(open('../DEC5D_decagon_order5_depth/step_53.py').read()[:3000])
