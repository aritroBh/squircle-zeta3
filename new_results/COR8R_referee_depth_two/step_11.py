import os, re
base = '<project>'
tex = open(os.path.join(base,'paper/ieee/main.tex')).read()
# find prop:M3angle
i = tex.find('prop:M3angle')
print(tex[max(0,i-3000):i+5000])
