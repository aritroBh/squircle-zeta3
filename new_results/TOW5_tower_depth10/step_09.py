with open('<project>/research/08_theory/TOW4_tower_depth8/tow2lib.py') as f:
    txt = f.read()
# find the truncated part
idx1 = txt.find('def complete_bell')
idx2 = txt.find('def tau1(n):')
print(txt[idx1:idx2])
