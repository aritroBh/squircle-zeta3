base = '<project>/research/08_theory/'
t = open(base+'TOW1_tower_law/proof_tower.md').read()
i2 = t.find('## 2'); i4 = t.find('## 4')
print(t[i2:i4])
