
import os
BASE = '<project>'
tex = open(os.path.join(BASE,'paper/squircle_note/main.tex')).read()
idx = tex.find('rem:k4')
print("=== around rem:k4 ===")
print(tex[idx-1500:idx+2500])
