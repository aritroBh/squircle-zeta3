print(open('../WTS_weight_k_structure/theory/run_kernel_k.py').read())
t = open('../WTS_weight_k_structure/theory/fastred.py').read()
import re
i = t.find('real')
print(t[max(0,i-2000):i+2000])
