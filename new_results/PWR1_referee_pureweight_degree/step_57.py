import os
# find D11c folder (operator representation K,T machinery)
for d in os.listdir('..'):
    if 'D11c' in d: print(d, os.listdir(f'../{d}')[:30])
