base = '<project>/research/08_theory'
k1e = f'{base}/K1e_hexagon_closed'
# read step_11 (6-arc check) to see how logs were computed from geometry.json
src = open(f'{k1e}/step_11.py').read()
print(src[:5000])
