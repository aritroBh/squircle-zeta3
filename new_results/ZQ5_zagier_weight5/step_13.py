import subprocess, sys
r = subprocess.run(['grep','-n','def emb\\|class QD\\|import mpmath\\|^import\\|^from\\|math.sqrt\\|def emb','w5lib.py'],capture_output=True,text=True)
print(r.stdout)
r2 = subprocess.run(['grep','-n','float\\|\\.sqrt','w5lib.py'],capture_output=True,text=True)
print(r2.stdout)
r3 = subprocess.run(['wc','-l','w5lib.py'],capture_output=True,text=True)
print(r3.stdout)
