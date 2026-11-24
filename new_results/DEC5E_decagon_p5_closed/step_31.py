import sys, pickle
sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 130
import gengine_dec2p
G1, letters = gengine_dec2p.gtab(120, 380, cachedir=".")
raw0, letters0 = pickle.load(open('../DEC5D_decagon_order5_depth/gtab_dec_d120_N380.pkl','rb'))
G0 = {k: mp.mpf(v) for k, v in raw0.items()}
worst = mp.mpf(0); wk=None
for w in G0:
    d = abs(G1[w]-G0[w])
    if d > worst: worst, wk = d, w
print("words:", len(G0), "worst diff vs DEC5D validated table:", mp.nstr(worst,3), "at", wk)
open('progress.log','a').write("step_29: gengine_dec2p validated vs DEC5D gtab_dec_d120_N380.pkl: worst diff %s over %d words\n" % (mp.nstr(worst,3), len(G0)))
