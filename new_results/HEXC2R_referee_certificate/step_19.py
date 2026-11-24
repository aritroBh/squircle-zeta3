
import subprocess, json
out = subprocess.run(["find", "<project>/research", "-name", "hpl_series.py"], capture_output=True, text=True)
print(out.stdout)
