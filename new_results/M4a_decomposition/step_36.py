import os
BASE = '<project>')
print(os.listdir(os.path.join(BASE,'paper/squircle_note/scripts')))
src = open(os.path.join(BASE,'paper/squircle_note/scripts/m3_tails_exact.py')).read()
print(src)
