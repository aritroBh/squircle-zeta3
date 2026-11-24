
src = open('step_37.py').read()
# print the truncated region: find "G1B  = a" up to "p2m = arb"
i = src.find('G1B  = a')
j = src.find('p2m = arb')
print(src[i-200:j])
