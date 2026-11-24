txt = open('../DCY2_decay_first_correction/dcy2_cell.py').read()
i = txt.find('def Gsplit_of')
print(txt[i:])
