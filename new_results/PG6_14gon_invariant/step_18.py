import os
base = '<project>/research/08_theory/'
dods = [d for d in os.listdir(base) if d.startswith('DOD')]
for d in sorted(dods): print(d)
print("---- searching for P3 / P''' known values in DOD result.md files ----")
import re
for d in sorted(dods):
    rm = os.path.join(base,d,'result.md')
    if os.path.exists(rm):
        txt = open(rm).read()
        for mline in txt.splitlines():
            if ("P'''" in mline or 'P3' in mline or "P‴" in mline) and ('-' in mline):
                print(f"[{d}] {mline.strip()[:150]}")
