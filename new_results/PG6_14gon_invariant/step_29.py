import re
txt = open('<project>/paper/ieee/main.tex').read()
for m in re.finditer(r'Polygons at the Diamond End', txt):
    print("SECTION at", m.start())
print("="*70)
print(txt[214900:216500])
