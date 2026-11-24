base = '<project>/research/08_theory/'
t = open(base+'TOW1_tower_law/proof_tower.md').read()
i5 = t.find('## 5')
i6 = t.find('## 6')
print(t[i5:i6])
