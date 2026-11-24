open('progress.log','a').write('step_22..29: W3 affine law EXACT (4 angles k<=14); slope kernels universal x s(alpha)=(1-c^2)/(2c^2); slope value = s*J(r), J resists classical basis -> computing exact HPL form\n')
import os
root = '<project>'
a04dir = root + '/research/10_frontier/W1_corner_family/A_function_level/'
print([f for f in os.listdir(a04dir) if 'hpl' in f.lower()])
print(open(a04dir + 'hpl_tools.py').read()[:3000] if os.path.exists(a04dir+'hpl_tools.py') else 'no hpl_tools')
