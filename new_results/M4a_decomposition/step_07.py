import os
BASE = '<project>'
tex = open(os.path.join(BASE,'paper/squircle_note/main.tex')).read()
print(tex[14000:17000])
