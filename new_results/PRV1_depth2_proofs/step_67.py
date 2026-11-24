wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/result.md') as f:
    txt = f.read()
txt = txt.replace("`step_01.py`–`step_62.py`", "`step_01.py`–`step_66.py`")
with open(f'{wd}/result.md','w') as f:
    f.write(txt)

# Final sanity: print key segments of verify_full.out (the part 3c + part 4 + ALL CHECKS)
with open(f'{wd}/verify_full.out') as f:
    out = f.read()
i = out.find('PART 3c')
print(out[i:i+1900])
