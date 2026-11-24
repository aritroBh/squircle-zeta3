import os, subprocess

# First, find the research directory
for base in ['<project>', os.getcwd()]:
    p = os.path.join(base, 'research/08_theory')
    if os.path.isdir(p):
        print("FOUND:", p)
        print(sorted(os.listdir(p)))
        break
else:
    # search
    for root, dirs, files in os.walk('<home>/Downloads/Research', topdown=True):
        if 'TOW3_tower_depth6' in dirs:
            print("found at", root)
            break
