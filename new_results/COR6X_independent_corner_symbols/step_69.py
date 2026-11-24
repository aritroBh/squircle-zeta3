import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c","""
txt = open('<project>/research/08_theory/COR6_true_point_symbols/result.md').read()
print(txt[:8000])
"""], capture_output=True, text=True)
print(r.stdout[:8500])
