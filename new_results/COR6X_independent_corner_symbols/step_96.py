import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",
  "print(open('<project>/research/08_theory/COR6_true_point_symbols/cor6_ft.py').read()[:4000])"],
  capture_output=True, text=True)
print(r.stdout[:4200])
