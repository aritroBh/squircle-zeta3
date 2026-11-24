
base='research/08_theory/DEC1_decagon'
txt = open(base+'/result.md').read()
i = txt.find('## 3')
j = txt.find('## 4')
print(txt[i:j])
