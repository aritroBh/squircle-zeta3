base = '<project>/research/08_theory/'
with open(base+'PG5_alphabet_14gon/result.md') as f:
    content = f.read()
i = content.find('s = m')
j = content.find('Task 4')
print(content[i-100:j])
