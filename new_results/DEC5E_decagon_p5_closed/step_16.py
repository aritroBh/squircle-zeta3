import subprocess
out = subprocess.run(['grep','-n','Li3\|Li_3\|lewin\|Lewin\|golden\|reduc','../DEC2_decagon_exact_assembly/verify_dec2.py'],capture_output=True,text=True)
print(out.stdout[:5000])