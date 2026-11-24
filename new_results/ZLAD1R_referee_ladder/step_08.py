import subprocess, os
base = '<project>'
print(os.path.exists(base + '/paper/ieee/main.tex'))
out = subprocess.run(['grep', '-n', 'thm:flat', base + '/paper/ieee/main.tex'], capture_output=True, text=True)
print(out.stdout)
