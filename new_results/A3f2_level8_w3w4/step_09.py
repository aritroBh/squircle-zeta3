import os
d='<project>/research/08_theory'
print(open(d+'/A3f2_level8_w3w4.brief').read()[:500])
print('=== dir A3f2_level8_w3w4:', os.listdir(d+'/A3f2_level8_w3w4'))
print('=== .run:', os.listdir(d+'/A3f2_level8_w3w4.run'))
print('=== A3f .run:', os.listdir(d+'/A3f_level8_relations.run'))
