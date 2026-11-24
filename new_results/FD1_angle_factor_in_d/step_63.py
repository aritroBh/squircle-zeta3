t2='<project>/research/08_theory/TOW2_tower_all_depths'
src=open(t2+'/tow2lib.py').read()
i=src.find('def rA0')
print(src[i:i+4700])
