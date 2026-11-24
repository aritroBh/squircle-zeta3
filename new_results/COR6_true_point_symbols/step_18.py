import subprocess
r = subprocess.run(["grep", "-n", "eq:R}", "<project>/paper/ieee/main.tex"], capture_output=True, text=True)
print(r.stdout)
r = subprocess.run(["grep", "-n", "label{eq:R}", "<project>/paper/ieee/main.tex"], capture_output=True, text=True)
print(r.stdout)
