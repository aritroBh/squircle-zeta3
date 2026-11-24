import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
import mpmath as mp
t0=time.time()
OCT5D = "<project>/research/08_theory/OCT5D_octagon_order5_depth/"
parsed = pickle.load(open(OCT5D+"parsed_mons.pkl","rb"))
print("parsed loaded", time.time()-t0)
rewrite = pickle.load(open("rewrite_all.pkl","rb"))
print("rewrite loaded", time.time()-t0, os.path.getsize("rewrite_all.pkl")/1e6, "MB")
import gengine2p
G1, letters = gengine2p.gtab(430, 1440, cachedir="<project>/research/08_theory/OCT5D2_octagon_depth_test")
print("gtab loaded", time.time()-t0)
