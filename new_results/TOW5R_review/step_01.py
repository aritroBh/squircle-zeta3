import os, subprocess
# First, find the research directory
for root in ['<project>', '.']:
    try:
        print(root, os.listdir(root)[:50])
    except Exception as e:
        print(root, 'ERR', e)
