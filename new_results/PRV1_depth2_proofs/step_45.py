wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/verify_full.py') as f:
    content = f.read()
i = content.find('def part0')
j = content.find('# Part 1')
print(content[i:j])
