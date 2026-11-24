
import os
BASE = '<project>'
# 1) paper main.tex: find thm:M3 and rem:k4
tex = open(os.path.join(BASE,'paper/squircle_note/main.tex')).read()
import re
for tag in ['thm:M3','rem:k4','conj:wt2','label{']:
    pass
idx = tex.find('thm:M3')
print("=== around thm:M3 ===")
print(tex[idx-3000:idx+3000])
