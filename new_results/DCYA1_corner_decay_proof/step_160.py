import sys, math
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
from arb_ei import E_I_sc, iv
def thin_cell(N, dw=0.1):
    s=math.sqrt(N); s2=s+dw
    n1=int(math.floor(s*s)); n2=int(math.ceil(s2*s2))
    if n2<=n1: n2=n1+1
    return iv(str(n1),str(n2))
for N in [100,200,400,800,1600,2400,3200,4800,6400,9600,12800,16000,20000,24000,32000,48000,64000,96000,128000,192000,256000,384000,512000,768000,1000000]:
    try:
        v=float(E_I_sc(thin_cell(N)))
        print(f"N={N}: E_I/sc={v:.2f}")
    except Exception as e:
        print(f"N={N}: ERR {type(e).__name__}")
