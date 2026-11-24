
base = "<project>/research/08_theory"
import subprocess
out = subprocess.run(["grep", "-l", "F4_relation", base + "/COR5_corner_classicality/" + f"step_{i}.py" for i in range(70, 120)], capture_output=True, text=True)
print(out.stdout, out.stderr)
