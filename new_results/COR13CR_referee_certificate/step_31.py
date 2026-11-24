import os
print(open('step_01.py').read())
print("---")
base = '<project>/research/08_theory'
print(open(os.path.join(base,'COR13C_third_corner_certificate','progress.log')).read())
print("---")
print(os.listdir('.'))
