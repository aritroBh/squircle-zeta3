src = open('../CHN2_chain_planes/mzv_exact.py').read()
# find BASIS section
i = src.find('BASIS = {')
print(src[i:i+1400])
