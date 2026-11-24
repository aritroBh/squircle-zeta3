import subprocess
out = subprocess.run(["grep","-n","12-gon delta","../PGEN1_polygon_alphabet/step_46.py","../PGEN1_polygon_alphabet/step_57.py","../PGEN1_polygon_alphabet/step_54.py"],capture_output=True,text=True).stdout
print(out)
src = open('../PGEN1_polygon_alphabet/step_46.py').read()
i = src.find('12-gon delta')
print(src[max(0,i-3500):i+500])
