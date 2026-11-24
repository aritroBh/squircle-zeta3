base = '<project>'
with open(base + '/paper/ieee/main.tex') as f:
    lines = f.readlines()
print(''.join(lines[1400:1470]))
