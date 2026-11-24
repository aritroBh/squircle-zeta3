with open('../../../paper/ieee/main.tex') as f:
    tex = f.read()
i = tex.find('Polygons at the Diamond End')
print(tex[i+9000:i+20000])
