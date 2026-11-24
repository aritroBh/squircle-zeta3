
import re
tex = open('<project>/paper/ieee/main.tex').read()
print("LENGTH:", len(tex))
# Find the relevant sections
for pat in [r'How Far Do the Series Reach', r'Radius Zero', r'Radius Exactly One', r'Superellipsoid', r'lem:holo', r'thm:oct', r'thm:S']:
    for m in re.finditer(pat, tex):
        print(pat, '->', m.start())
