import os, time
work = '<project>/research/08_theory/TOW4_tower_depth8'
log = open(work + '/progress.log','a')
def logline(s):
    log.write(time.strftime('%H:%M:%S ') + s + '\n'); log.flush()
    print(s)
logline("SETUP: copied tow2lib/tow3lib + L/arho/e61/e_fits caches + prove_*.py from TOW3 (byte copies).")
logline("TASK A: a_7 = N_7/(40320 Pi_3), a_8 = N_8/(645120 Pi_4) assembled via Bell from L(1..8) cache (2.0s, 18.1s).")
logline("TASK A: rho(7)/e00==a_7, rho(8)/e00==a_8 vs TOWN1 towers (k0=d all keys, no mult): 717 cells each, 0 bad.")
logline("TASK A: L_simplified cache == tow2lib.L(r) proved construction SYMBOLICALLY, r=1..8 (diff==0; r=8 9.9s). Closes TOW3's r<=6 gap at r=7,8.")
log.close()
print("logged")
