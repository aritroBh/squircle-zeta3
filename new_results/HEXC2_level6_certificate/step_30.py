src2 = open('hexc2_common.py').read()
i = src2.find('def li_atom')
j = src2.find('f1 =', i)
print(src2[j-120:j+220])
