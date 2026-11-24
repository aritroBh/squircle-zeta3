import os, subprocess
base = 'research/08_theory'
print(os.getcwd())
print(sorted(os.listdir(base)))
# create work dir
work = os.path.join(base, 'TOW1_tower_law')
os.makedirs(work, exist_ok=True)
print('work dir:', work, os.path.isdir(work))
