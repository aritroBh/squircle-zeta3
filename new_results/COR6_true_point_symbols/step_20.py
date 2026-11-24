import subprocess
r = subprocess.run(["grep", "-n", "Phi(c)\\|Psi(c)\\|Phi\\|Psi", "<project>/paper/ieee/main.tex"], capture_output=True, text=True)
print(r.stdout[:4000])
