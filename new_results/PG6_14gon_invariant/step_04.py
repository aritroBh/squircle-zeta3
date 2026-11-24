with open(base+'PG5_alphabet_14gon/result.md') as f:
    content = f.read()
# Print the section between "s = m−2 = 5" and "δ(N₃) verdict"
i = content.find('s = m')
j = content.find('δ(N₃) verdict')
print(content[i-200:j+200])
