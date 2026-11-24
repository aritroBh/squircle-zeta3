
import os
BASE = '<project>'
tex = open(os.path.join(BASE,'paper/squircle_note/main.tex')).read()
idx = tex.find('label{rem:k4}')
print(tex[idx-1200:idx+2200])
