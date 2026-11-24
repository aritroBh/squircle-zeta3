import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
print(open("build_big.log").read())
print("cache files:", [f for f in os.listdir(".") if f.startswith("gtab2_")])
