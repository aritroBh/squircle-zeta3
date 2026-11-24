import subprocess
out = subprocess.run(["grep","-rln","L3d\\|3397\\|oct_delta\\|octagon delta","--include=step_*.py","../PGEN1_polygon_alphabet"],capture_output=True,text=True).stdout
print(out)
