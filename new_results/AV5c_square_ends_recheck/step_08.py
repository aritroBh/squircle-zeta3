txt = open('../AV5_layer_universality/result.md').read()
i = txt.find('Octagon')
print(txt[i:i+3000])
