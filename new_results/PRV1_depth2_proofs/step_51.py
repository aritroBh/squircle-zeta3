p = 'research/08_theory/DEC1_decagon/result.md'
with open(p) as f:
    txt = f.read()
i = txt.find('§3') if '§3' in txt else txt.find('## 3')
if i < 0:
    i = txt.find('Li_{2,1}')
print(txt[max(0,i-200): i+3000])
