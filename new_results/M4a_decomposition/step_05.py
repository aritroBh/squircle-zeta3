import os, datetime
BASE = '<project>'
with open(os.path.join(BASE,'paper/squircle_note/main.tex')) as f:
    tex = f.read()
# Print the area around prop:split and W4 definition
print(tex[5600:10000])
