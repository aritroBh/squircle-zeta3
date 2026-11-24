import os, subprocess, textwrap, sys, json, pathlib

BASE = '<project>'
# Explore the structure
for path in ['paper/squircle_note', 'paper/squircle_note/scripts', 'paper/squircle_note/outputs']:
    full = os.path.join(BASE, path)
    print('=== ', path)
    try:
        for f in sorted(os.listdir(full)):
            fp = os.path.join(full, f)
            print(f, os.path.getsize(fp))
    except Exception as e:
        print('ERR', e)
