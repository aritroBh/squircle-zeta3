import sys, pickle
sys.path.insert(0,'.')
import mpmath as mp
# verify Li32 convention against G-table at 800
mp.mp.dps = 800
import gengine_dec2p
G1, letters = gengine_dec2p.gtab(800, 2700, cachedir=".")
L0,L1,LP,LQ = letters['L0'],letters['L1'],letters['LP'],letters['LQ']
li32v = pickle.load(open('li32_gold_800.pkl','rb'))
phi = (1+mp.sqrt(5))/2
v = mp.mpf(li32v['vals'][0])
print("Li32(-1/phi,1) vs G(0,0,-phi,0,-phi;1):", mp.nstr(abs(v - G1[(L0,L0,LP,L0,LP)]), 3))
v1 = mp.mpf(li32v['vals'][1])
print("Li32(1/phi^2,1) vs G(0,0,phi^2,0,phi^2;1):", mp.nstr(abs(v1 - G1[(L0,L0,LQ,L0,LQ)]), 3))
v2 = mp.mpf(li32v['vals'][2])
print("Li32(-1/phi,-1/phi) vs G(0,0,-phi,0,phi^2;1):", mp.nstr(abs(v2 - G1[(L0,L0,LP,L0,LQ)]), 3))
v3 = mp.mpf(li32v['vals'][3])
print("Li32(1/phi^2,-phi) vs G(0,0,phi^2,0,-phi;1):", mp.nstr(abs(v3 - G1[(L0,L0,LQ,L0,LP)]), 3))
