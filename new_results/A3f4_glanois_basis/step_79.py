import subprocess
base = '<home>)'
for pat in ['22-term','22 term','twenty-two','Goncharov','22_term']:
    out = subprocess.run(['grep','-rli',pat, base+'/research'],capture_output=True,text=True).stdout
    print(pat, '->', out)
