import os, re
paper = '<project>/paper/ieee/main.tex'
with open(paper) as f:
    txt = f.read()
# Find prop:M3angle
idx = txt.find('prop:M3angle')
while idx != -1:
    print('--- occurrence at', idx, '---')
    print(txt[max(0,idx-200):idx+3000])
    print()
    idx = txt.find('prop:M3angle', idx+1)
    if idx > 20000: break
