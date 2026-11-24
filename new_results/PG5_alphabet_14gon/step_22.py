import subprocess, os, glob
# find DOD1, DOD2, DOD3 directories
for d in sorted(os.listdir('..')):
    if 'DOD' in d or 'dod' in d.lower():
        print(d)
print()
src = open('../../08_theory/D17W_draft17_writer/mirror_backup/main.tex').read()
# the 12-gon paragraph
i = src.find('12-gon')
while i != -1 and i < 145000:
    i = src.find('12-gon', i+1)
    if 139000 < i < 146000: break
print(src[139600:i+3000] if i!=-1 else src[139600:144000])
