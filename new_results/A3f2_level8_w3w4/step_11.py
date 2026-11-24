import os
d='<project>/research/08_theory'
print('cwd:', os.getcwd())
print(sorted(os.listdir('.')))
print('A3f .run is file:', os.path.isfile(d+'/A3f_level8_relations.run'))
