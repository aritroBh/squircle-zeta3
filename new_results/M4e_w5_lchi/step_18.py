import os
p = '<project>/paper/squircle_note/scripts'
print(sorted(os.listdir(p)))
src = open(os.path.join(p,'qsqrt2.py')).read()
print(src[:4000])
