src = open('gold5_lib.py').read()
i = src.find('def quotient_coords_nf')
print(src[i:i+900])
