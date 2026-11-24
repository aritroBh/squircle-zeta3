with open('<project>/research/08_theory/TOW2_tower_all_depths/proof_all_depths.md') as f:
    txt = f.read()
i1 = txt.find('## 1')
i2 = txt.find('## 2')
print(txt[i1:i2])
