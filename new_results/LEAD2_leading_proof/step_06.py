lead1 = '<project>/research/08_theory/LEAD1_leading_genfun'
with open(lead1 + '/genfun.md') as f:
    txt = f.read()
# print the middle portion that got truncated (sections 3.2 and 3.3)
i = txt.find('### 3.2')
j = txt.find('Hence for even d')
print(txt[i:j])
