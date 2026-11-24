import subprocess
r=subprocess.run(['<project>/.venv/bin/python','gen_supplement.py'],capture_output=True,text=True)
print(r.stderr[-800:])
g=open('gen_supplement.py').read()
i=g.find("A(r'Numerically')")
print(g[i-200:i+900])
