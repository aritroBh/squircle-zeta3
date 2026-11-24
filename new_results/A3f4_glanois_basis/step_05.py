p = 'research/08_theory/A3g_denominator/glanois_1411.4947.txt'
with open(p) as f:
    lines = f.readlines()
print(len(lines))
print(''.join(lines[180:260]))
