import os
# find the step files that mention key results
for f in sorted(os.listdir('oct5')):
    if f.startswith('step'):
        src = open(os.path.join('oct5',f)).read()
        if 'exact_LS' in src or 'mixed' in src or 'atoms' in src or 'Jloc' in src or 'u =' in src or 'P3_720' in src:
            print('='*30, f, '='*30)
            print(src[:2500])
