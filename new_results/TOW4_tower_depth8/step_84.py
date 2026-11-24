import time
work = '<project>/research/08_theory/TOW4_tower_depth8'
for _ in range(90):
    txt = open(work + '/verify_tow4.log').read()
    if 'VERIFY_TOW4' in txt: break
    time.sleep(6)
print(txt)
