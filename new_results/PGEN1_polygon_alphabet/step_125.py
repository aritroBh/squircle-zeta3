# step_96: capture full verify output
import subprocess
r = subprocess.run(["nice", "-n", "10", "<project>/.venv/bin/python", "verify_pgen1.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/PGEN1_polygon_alphabet")
open("verify_output.txt","w").write(r.stdout)
i0 = r.stdout.find("PART 1")
i1 = r.stdout.find("PART 2")
print(r.stdout[i0:i1])
