t = r'''
import re, json
txt = open("<project>/paper/ieee/main.tex").read()
i = txt.find("sec:corners-B")
# search for the M3(2pi/3) closed form near "Li_4" and "1/9"
for m in re.finditer(r"frac19\}|tfrac19\}|\{1/9\}", txt):
    s = txt[m.start()-1500:m.start()+800]
    if "Li_4" in s:
        print(s)
        print("=====================")
# Psi file
psi = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(1,0,0).json"))
print("Psi file:", len(psi), "terms")
for s, c in psi:
    print(s, "|", c)
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t], capture_output=True, text=True)
print(r.stdout[:9000])
