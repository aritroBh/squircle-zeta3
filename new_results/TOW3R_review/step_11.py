base = '<project>/research/08_theory/'
txt = open(base+'CX8_tow2_review/report.md').read()
i = txt.find('## 5')
print(txt[i:i+7000])
