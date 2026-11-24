import os, re
# find the paper
for cand in ['../../paper/squircle_note/main.tex', '../../../paper/squircle_note/main.tex']:
    print(cand, os.path.exists(cand))
# locate TOW2
base = '..'
for f in os.listdir(base):
    if 'TOW' in f.upper():
        print(f)
