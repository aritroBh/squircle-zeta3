import subprocess
t=open('supplement_W5.tex').read()
i=t.find('classical part')
print(t[i-100:i+700])
r=subprocess.run(['nice','-n','10','/opt/homebrew/bin/tectonic','-X','compile','supplement_W5.tex'],capture_output=True,text=True,timeout=600)
print('--- tectonic supplement ---'); print(r.stdout[-1500:]); print(r.stderr[-1500:])
