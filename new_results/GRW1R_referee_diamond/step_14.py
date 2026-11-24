import json
d = json.load(open('../D11b2_stokes_jump/step2_num450_N180.json'))
f = d['f450']
for i in range(8): print(i, f[i][:50])
# find how GRW1 reads it
src = open('../GRW1_diamond_growth_borel/step_04.py').read()
print(src[:3000])
