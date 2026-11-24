import os
wd = '<project>/research/08_theory/TOW5R_review'
os.chdir(wd)
print(open('progress.log').read())
print('---files---')
print(sorted(os.listdir('.')))
