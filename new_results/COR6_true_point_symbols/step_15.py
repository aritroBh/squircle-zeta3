base = "<project>/research/"
import subprocess, glob
# find paper
for pat in ["paper/ieee/main.tex", "../paper/ieee/main.tex"]:
    p = base + pat
    r = subprocess.run(["ls", p], capture_output=True, text=True)
    print(p, r.stdout, r.stderr)
r = subprocess.run(["find", "<project>", "-name", "main.tex", "-path", "*ieee*"], capture_output=True, text=True)
print(r.stdout)
