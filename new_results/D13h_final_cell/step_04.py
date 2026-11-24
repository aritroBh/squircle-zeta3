import re
p = 'research/08_theory/D13g_all_N/result.md'
txt = open(p).read()
print("total chars:", len(txt))
# Find section 4
i = txt.find('## 4')
if i < 0: i = txt.find('# 4')
if i < 0: i = txt.find('§4')
# print from there to next ## 5 or end
j = txt.find('## 5', i+1)
if j < 0: j = min(i+12000, len(txt))
print(txt[i:j][:12000])
