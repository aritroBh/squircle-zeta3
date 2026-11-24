import os
# set up work dir
wd='research/08_theory/D13t_adversarial_review'
os.makedirs(wd, exist_ok=True)
os.chdir(wd)
open('progress.log','w').write('# D13t adversarial review progress log\n')
print(os.getcwd())
# find D13r review
base='research/08_theory/'
for root,dirs,files in os.walk(base):
    if 'D13r' in root:
        for f in sorted(files):
            print(os.path.join(root,f))
