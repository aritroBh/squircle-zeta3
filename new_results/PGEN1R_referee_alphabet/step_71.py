src = open('../PGEN1_polygon_alphabet/step_46.py').read()
i = src.find('12-gon delta')
# find the data table for the 12-gon coefficients
import re
m = re.search(r'data12\s*=\s*', src)
print(src[m.start()-200:m.start()+1200] if m else 'no data12')
