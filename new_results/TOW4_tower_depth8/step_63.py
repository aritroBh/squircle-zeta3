import os, time, pickle
base = '<project>/research/08_theory'
work = base + '/TOW4_tower_depth8'
log = open(work + '/progress.log','a')
log.write(time.strftime('%H:%M:%S ') + "FULL LAW delta<=8 vs TOWN1: 49953 cells (1596 diagonal), 0 bad; +222 diagonal-initial cells, 0 bad. Tower law PROVED to depth 8.\n")
log.close()
wts = base + '/WTS_weight_k_structure'
print("w6:", sorted(os.listdir(wts + '/w6')))
print("w5:", sorted(os.listdir(wts + '/w5'))[:20])
