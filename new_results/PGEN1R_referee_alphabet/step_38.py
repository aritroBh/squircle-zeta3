import subprocess
out = subprocess.run(["grep","-rln","shuffle","--include=step_*.py","../PGEN1_polygon_alphabet"],capture_output=True,text=True).stdout
print(out)
out2 = subprocess.run(["grep","-rn","shuffle","--include=*.md","../PGEN1_polygon_alphabet"],capture_output=True,text=True).stdout
print(out2[:2000])
