txt = open('../OCT9b_octagon_debug/result_nofinal.md').read()
# print sections 4 and 5
i4 = txt.find('## 4')
i6 = txt.find('## 6')
print(txt[i4:i6])
