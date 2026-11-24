
with open('step_03.py') as f:
    txt = f.read()
txt = txt.replace("A1=np.array([0,0,pA]); P1=np.array([pB,pB,pB]); Q1=np.array([pB,-pB,pB]); R1=np.array([-pB,-pB,pB])",
                  "A1=np.array([0.,0.,float(pA)]); P1=np.array([float(pB)]*3); Q1=np.array([float(pB),-float(pB),float(pB)]); R1=np.array([-float(pB),-float(pB),float(pB)])")
txt = txt.replace("A2=np.array([pA,0,0])", "A2=np.array([float(pA),0.,0.])")
with open('step_03.py','w') as f:
    f.write(txt)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','step_03.py'],
                   capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-3000:] if r.returncode else "exit 0")
