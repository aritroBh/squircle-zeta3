with open('<project>/research/08_theory/TOW4_tower_depth8/result.md') as f:
    txt = f.read()
idx = txt.find('## Task D')
idx2 = txt.find('## Evidence table')
print(txt[idx:idx2])
