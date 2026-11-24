import subprocess
r = subprocess.run(["grep", "-n", "thm:W4R", "<project>/paper/ieee/main.tex"], capture_output=True, text=True)
print(r.stdout)
