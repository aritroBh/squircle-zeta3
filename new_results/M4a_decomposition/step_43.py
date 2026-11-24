import os
BASE = '<project>'
OUT = os.path.join(BASE, 'paper/squircle_note/outputs')
print(open(os.path.join(OUT,'m3_closed_form.out')).read())
print('='*40)
print(open(os.path.join(OUT,'m3_assemble.out')).read())
