import pickle, sys
from fractions import Fraction as Fr
from math import comb
sys.path.insert(0,'.')
import ftfield as F
E3 = pickle.load(open('E3n14.pkl','rb'))['E3n']
constvec = pickle.load(open('constvec14.pkl','rb'))
# ORDER of term-logs in E3n keys: [3,4,1,2,5,6,0]
# term -> (lettervec, wlog_index) where wlog_index: 'lt','l1t','w1','w2','w5','w6', or None (const)
TERMS = ['a3','a4','a1','a2','a5','a6','a0']
WL = {'a3':'lt','a4':'l1t','a1':'w1','a2':'w2','a5':'w5','a6':'w6','a0':None}
ZEROV=(0,0,0)
def vadd(u,v): return (u[0]+v[0],u[1]+v[1],u[2]+v[2])
def vscale(v,k): return (v[0]*k,v[1]*k,v[2]*k)
# mons accumulator: key (letters, a, b, wc4) -> dict n->FT
mons = {}
def add_term(letters, a, b, wc, n, coefFT):
    key = (letters, a, b, wc)
    if key not in mons: mons[key]={}
    mons[key][n] = F.fadd(mons[key].get(n,F.ZERO), coefFT)
for key7, poly in E3.items():
    # start with (letters=ZEROV, wlogs dict, t-poly dict {n:FT})
    # expand product over 7 terms
    state = [ (ZEROV, {'lt':0,'l1t':0,'w1':0,'w2':0,'w5':0,'w6':0}) ]
    for ti,term in enumerate(TERMS):
        e = key7[ti]
        if e==0: continue
        lv = constvec[term]; wl = WL[term]
        newstate = []
        for (Lvec, W) in state:
            if wl is None:
                # only letter power, no binomial
                newstate.append( (vadd(Lvec, vscale(lv,e)), dict(W)) )
            else:
                for k in range(e+1):
                    W2 = dict(W); W2[wl]+=k
                    newstate.append( (vadd(Lvec, vscale(lv,e-k)), W2, comb(e,k)) )
        if wl is None:
            state = [(L,W) for (L,W) in newstate]
        else:
            state = [(L,W,c) for (L,W,c) in newstate]
        # normalize: carry binomial coeff
        if ti==0 or True:
            pass
    # state now list of (Lvec, W) or (Lvec,W,c)... handle uniformly below
    # redo cleanly:
print("approach too messy; redoing cleanly")
