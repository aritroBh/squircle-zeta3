import json, mpmath as mp
mp.mp.dps = 100
m5b='<project>/research/08_theory/M5b_w6_closed_form'
p46=json.load(open(m5b+'/w5pp_4600_partial.json'))
p48=json.load(open(m5b+'/w5pp_4800_partial.json'))
for k in ['left','mid','right']:
    x=mp.mpf(p46[k][0]); y=mp.mpf(p48[k][0])
    print(k, 'diff', mp.nstr(abs(x-y),5), mp.nstr(x,25))
# independent mpmath piece values
import json as js
om=js.load(open(m5b+'/omega5.json'))
r2=mp.sqrt(2); aN,bN=1+r2,r2-1
pole={"0":mp.mpf(0),"a":aN,"-a":-aN,"b":bN,"-b":-bN}
def piecequad(i,j,lo,hi):
    ts=om[f"{i},{j}"]
    def f(v):
        l1=mp.log((v-bN)*(v+aN)/(4*v)); l2=mp.log((aN-v)*(v+bN)/(4*v))
        return l1**i*l2**j*mp.fsum((mp.mpf(c[0])+r2*mp.mpf(c[1]))/(v-pole[z]) for z,c in ts)
    return mp.quad(f,[lo,(lo+hi)/2,hi])
for name,lo,hi in [('left',bN,1),('mid',1,aN-1),('right',aN-1,aN)]:
    s=mp.fsum([piecequad(*ij,lo,hi) for ij in [(2,2),(3,1),(4,0)]])
    print('mpmath',name,mp.nstr(s,25))
