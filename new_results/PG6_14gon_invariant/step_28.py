txt = open('<project>/paper/ieee/main.tex').read()
# Find the Polygons at the Diamond End section
for m in re.finditer(r'Polygons at the Diamond End', txt):
    print("SECTION at", m.start())
# Print around the (2,1) invariant hits
print("="*70)
print(txt[214900:216400])
