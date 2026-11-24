root = '<project>'
src = open(root + '/research/10_frontier/W1_corner_family/A_function_level/a04_rfamily.py').read()
i = src.find('Uc = to_hpl')
print(src[i:i+3500])
