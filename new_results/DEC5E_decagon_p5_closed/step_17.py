src = open('../DEC2_decagon_exact_assembly/verify_dec2.py').read()
i = src.find('PART 6')
print(src[i:i+7000])