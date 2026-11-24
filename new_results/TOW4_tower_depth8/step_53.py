import time
work = '<project>/research/08_theory/TOW4_tower_depth8'
for _ in range(60):
    txt = open(work + '/adj78_all.log').read()
    if 'ALL:' in txt: break
    time.sleep(5)
print(txt)
