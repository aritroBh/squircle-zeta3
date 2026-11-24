import subprocess
print(subprocess.run(["grep","-rn","shuffle","../PGEN1_polygon_alphabet/verify_pgen1.py"],capture_output=True,text=True).stdout[:3000])
