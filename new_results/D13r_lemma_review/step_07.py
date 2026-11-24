base='research/08_theory'
import os
print(os.listdir(base+'/D13e_certify'))
print('=====')
t=open(base+'/D13e_certify/result.md').read()
print(t[:5500])
