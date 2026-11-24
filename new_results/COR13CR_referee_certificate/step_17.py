import os, json
base = '<project>/research/08_theory'
c7 = os.path.join(base,'COR7_classical_corners')
print(open(os.path.join(c7,'c13_M3_quad_1100.txt')).read()[:800])
print()
w3 = json.load(open(os.path.join(c7,'c13_W3_relation.json')))
print("W3:", {k: v for k, v in w3.items() if k != 'coeffs' and k != 'basis'})
print("W3 a0:", w3['a0']); print("W3 basis:", w3['basis']); print("W3 coeffs:", w3['coeffs'])
f4 = json.load(open(os.path.join(c7,'c13_F4_relation.json')))
print("F4 a0:", f4['a0'])
print("F4 basis:", f4['basis'])
print("F4 coeffs:", f4['coeffs'])
