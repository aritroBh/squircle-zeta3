import subprocess
for i in [70,71,72]:
    src = open(f'../PGEN1_polygon_alphabet/step_{i:02d}.py').read()
    print('='*25, f'step_{i:02d}', '='*25)
    print(src[:2500])
