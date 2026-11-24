import os
BASE = '<project>'
tex = open(os.path.join(BASE,'paper/squircle_note/main.tex')).read()
print(tex[29800:31500])
print('===================== sec:M3 area =====================')
print(tex[32600:34800])
