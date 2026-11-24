import os
print('cwd:', os.getcwd())
print('parent (08_theory):', sorted(os.listdir('..')))
print('grandparent:', sorted(os.listdir('../..')))
