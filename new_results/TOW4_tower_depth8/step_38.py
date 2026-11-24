import time, os
work = '<project>/research/08_theory/TOW4_tower_depth8'
# poll until done (look for VERDICT) up to ~150s
for _ in range(30):
    txt = open(work + '/forced78.log').read()
    if 'VERDICT' in txt: break
    time.sleep(5)
print(txt[txt.find('delta=7 done'):])
