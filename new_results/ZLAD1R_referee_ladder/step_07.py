base = '<project>'
import subprocess
out = subprocess.run(['grep', '-n', 'thm:flat', '-r', base + '/paper'], capture_output=True, text=True)
print(out.stdout[:3000])
