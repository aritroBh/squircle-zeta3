work = '<project>/research/08_theory/M5a_next_coefficient'
print(open(f'{work}/progress.log').read())
print('='*80)
W6s = open(f'{work}/W6_400.txt').read().strip()
M5s = open(f'{work}/M5_400.txt').read().strip()
print("W6_400.txt (%d chars):" % len(W6s))
print(W6s)
print("M5_400.txt (%d chars):" % len(M5s))
print(M5s)