import os, re
base = "<project>/research/"
# search for the M4 value in M4e and M4a dirs
for d in ["08_theory/M4e_w5_lchi"]:
    p = base + d
    for f in sorted(os.listdir(p)):
        if f.endswith(('.json', '.md', '.txt')) and ('assemble' in f or 'closed' in f or 'result' in f):
            try:
                t = open(os.path.join(p, f)).read()
                for m in re.finditer(r'0\.0060279931\d{20,}', t):
                    print(f, ":", m.group(0)[:90])
                    break
            except Exception:
                pass
# also M4a
for d in os.listdir(base + "08_theory"):
    if d.startswith("M4a"):
        print("M4a dir:", d, sorted(os.listdir(base + "08_theory/" + d))[:30])
