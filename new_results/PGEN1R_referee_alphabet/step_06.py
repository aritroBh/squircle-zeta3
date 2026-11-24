import re
with open('../../../paper/ieee/main.tex') as f:
    tex = f.read()
# find the subsection
i = tex.find('Polygons at the Diamond End')
print(i, len(tex))
print(tex[i-200:i+9000])
