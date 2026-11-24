import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
import mpmath as mp, json
mp.mp.dps = 60
# paper R_1 formula
def P1_paper(m):
    V = 1/mp.tan(mp.pi/(2*m))
    P1 = 4*m*(1/mp.cos(mp.pi/(2*m)) - mp.cos(mp.pi/(2*m)))
    R1 = mp.fsum(mp.sin(i*mp.pi/m)*mp.log(V/mp.sin(i*mp.pi/m)) for i in range(1, m))/V
    return P1*R1
mine = {4: mp.mpf("2.75023713004966371870863005296385194252854498"),
        5: mp.mpf("2.72735150517983661148798551464963366249253"),
        6: mp.mpf("2.61219423683336630117516166502265676893478")}
for m in [4,5,6]:
    print(f"m={m}: P'(1) mine={mp.nstr(mine[m],30)}")
    print(f"      paper R_1 ={mp.nstr(P1_paper(m),30)}  resid={mp.nstr(mine[m]-P1_paper(m),4)}")
# m=4: P'(1) = u[(40-28sqrt2)log2 + (24-16sqrt2)log(1+sqrt2)]
u = mp.sqrt(2+mp.sqrt(2))
p1m4 = u*((40-28*mp.sqrt(2))*mp.log(2) + (24-16*mp.sqrt(2))*mp.log(1+mp.sqrt(2)))
print("m=4 paper explicit:", mp.nstr(p1m4, 30), " resid:", mp.nstr(mine[4]-p1m4, 4))
