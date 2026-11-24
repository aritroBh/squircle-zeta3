import os
work='research/08_theory/D13r_lemma_review'
print(open(work+'/step_01.py').read()[:800])
print('---2---')
print(open(work+'/step_02.py').read()[:800])
print([f for f in os.listdir(work)])
