with open('../../../paper/ieee/main.tex') as f:
    tex = f.read()
i = tex.find('Polygons at the Diamond End')
# find prop:edgelaw definition
j = tex.find('label{prop:edgelaw}')
print(tex[j-1500:j+2500])
