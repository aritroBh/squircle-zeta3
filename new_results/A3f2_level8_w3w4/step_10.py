import os
d='<project>/research/08_theory'
print('A3f .run:', os.listdir(d+'/A3f_level8_relations.run'))
print('A3f2 .run file:', open(d+'/A3f2_level8_w3w4.run').read()[:200] if os.path.isfile(d+'/A3f2_level8_w3w4.run') else 'is dir')
print('cwd:', os.getcwd())
print(os.listdir('.'))
