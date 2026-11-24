base = '<project>/research/08_theory'
t2 = f'{base}/TOW2_tower_all_depths'
src = open(f'{t2}/towlib.py').read()
print(len(src))
print(src[:6000])
