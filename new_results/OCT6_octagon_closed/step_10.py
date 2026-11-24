import os
for f in sorted(os.listdir('oct5')):
    if f.startswith('step'):
        src = open(os.path.join('oct5',f)).read()
        if 'rep' in src or 'Li3(' in src or 'Landen' in src or 'inversion' in src:
            print('='*28, f, '='*28)
            print(src)
            print()
