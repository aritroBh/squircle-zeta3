txt = open('../OCT9_generic_octagon_closed/result.md').read()
i = txt.find('Z_e')
print(txt[max(0,i-2500):i+3000])
