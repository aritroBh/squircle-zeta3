import os
base = '<project>/research'
for d in sorted(os.listdir(base)):
    if 'K2' in d or 'D13' in d:
        print(d)
print('---08_theory---')
base2 = os.path.join(base,'08_theory')
for d in sorted(os.listdir(base2)):
    if 'K2' in d or 'D13' in d:
        print(d)
