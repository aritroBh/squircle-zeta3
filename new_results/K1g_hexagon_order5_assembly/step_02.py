import os

base = '<project>'
k1f = os.path.join(base, 'research/08_theory/K1f_hexagon_order5')
print("=== Contents of 08_theory ===")
theory = os.path.join(base, 'research/08_theory')
for d in sorted(os.listdir(theory)):
    print(d)
print()
print("=== Contents of K1f ===")
for d in sorted(os.listdir(k1f)):
    print(d)
