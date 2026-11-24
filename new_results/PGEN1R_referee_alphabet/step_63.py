import subprocess
out = subprocess.run(["cat","../PGEN1_polygon_alphabet/verify_output.txt"],capture_output=True,text=True).stdout
print(len(out))
# print sections mentioning octagon/12-gon/decagon delta
import re
for m in re.finditer(r'(octagon|12-gon|decagon|delta)[^\n]*', out):
    pass
print(out[:5500])
