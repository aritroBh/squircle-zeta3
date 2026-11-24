import os
print(sorted(os.listdir('../COR6_true_point_symbols')))
src = open('../COR6_true_point_symbols/cor6_sunit.py').read()
print(src[:4000])
