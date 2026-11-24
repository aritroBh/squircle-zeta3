wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/verify_full.py') as f:
    content = f.read()
i = content.find('# Part 3: the decagon')
if i < 0:
    # find where decagon part starts
    i = content.find('G_dec')
    print(content[i-2000:i])
else:
    print(content[i:i+4200])
