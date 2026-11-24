with open('<project>/research/08_theory/TOW4_tower_depth8/proof_depth8.md') as f:
    txt = f.read()
# print lines 50-120
lines = txt.split('\n')
print("total lines:", len(lines))
print('\n'.join(lines[48:110]))
