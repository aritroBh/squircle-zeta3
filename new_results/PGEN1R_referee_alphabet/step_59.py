import subprocess
out = subprocess.run(["grep","-n","L3d\\|L3k\\|L3mr2h\\|L3r2h\\|L3x1\\|L3x2","../OCT6_octagon_closed/result.md"],capture_output=True,text=True).stdout
print(out[:3000])
