
s7 = open("steps/step_07.py").read()
s7 = s7.replace("from flint import fmpq_mat\nimport numpy as np",
                "from flint import fmpq_mat, fmpq\nimport numpy as np")
s7 = s7.replace("            M[j, len(rows)] = D.get(w, Fr(0))",
                "            M[j, len(rows)] = fmpq(str(D.get(w, Fr(0))))")
open("steps/step_07.py", "w").write(s7)
import subprocess
r = subprocess.run(["nice", "-n", "10", "<project>/.venv/bin/python", "steps/step_07.py"], capture_output=True, text=True, timeout=3500)
print(r.stdout[-1800:])
print(r.stderr[-3000:] if r.returncode else "OK")
