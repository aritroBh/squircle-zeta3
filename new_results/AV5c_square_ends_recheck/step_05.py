txt = open('../AV5_layer_universality/result.md').read()
print(len(txt))
# print section about square end
i = txt.find('Square end')
print(txt[i-200:i+4000])
