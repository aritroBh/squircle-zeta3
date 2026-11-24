import sys, os, pickle
sys.path.insert(0, os.getcwd())
from k_load import *
sys.path.insert(0, TH + '/TOWN1_tower_large_n')
import kernel_from_tower as KT
from kernel_from_tower import gm_from_tower, simple_part_of
from fastred import (QSI, ZERO, pf_add, pf_scale, pf_mul_DL, pf_fold, pf_split, pf_zero, DL1, DL2)
from fractions import Fraction as Fr
H = pickle.load(open(TH+'/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb')); H.pop('_meta')

def hdict(key, k, override=None):
    n=sum(key)
    out={}
    for dl in range(0, 2*n-3):
        t = (override or {}).get(dl, H[key][dl][k]) if override is not None else H[key][dl][k]
        out[2*n-2-dl]=QSI(*t)
    return out

def kernels_deg5(towers):
    """towers: dict key-> h dict {r:QSI} for keys (6,0),(5,1),(4,2),(3,3). returns W (full PF) for (5,0),(4,1),(3,2)"""
    g = {key: gm_from_tower(h) for key,h in towers.items()}
    W50 = pf_add(pf_scale(pf_mul_DL(g[(6,0)], DL1), Fr(-6)), pf_scale(pf_mul_DL(g[(5,1)], DL2), Fr(-1)))
    W41 = pf_add(pf_scale(pf_mul_DL(g[(5,1)], DL1), Fr(-5)), pf_scale(pf_mul_DL(g[(4,2)], DL2), Fr(-2)))
    W32 = pf_add(pf_add(pf_scale(pf_mul_DL(g[(4,2)], DL1), Fr(-4)),
                        pf_scale(pf_mul_DL(g[(3,3)], DL2), Fr(-3))),
                 pf_scale(pf_fold(pf_mul_DL(g[(3,3)], DL1)), Fr(-3)))
    return {(5,0):W50, (4,1):W41, (3,2):W32}

def realdict(P):
    return {(r.p0, r.p1): (c.p0,c.p1,c.p2,c.p3) for r,c in P[0].items()}
def stored_real(k,key):
    return {p:tuple(c) for p,c in kern(k)[key]['real'].items()}
def towers_at(k, over33=None):
    T={}
    for key in [(6,0),(5,1),(4,2),(3,3)]:
        T[key]=hdict(key,k)
    if over33 is not None:
        T[(3,3)]={10-dl:QSI(*over33[dl]) for dl in range(9)}
    return T

if __name__=="__main__":
    print("== Lemma-2 map rebuilt by me: kernel(deg5) from TOWN1 deg-6 towers vs stored full-reducer kernels")
    tot=0;bad=0
    for k in range(6,22):
        W = kernels_deg5(towers_at(k))
        for key in KEYS5:
            simple, higher, const = pf_split(W[key])
            mine = realdict(simple)
            mine = {p:c for p,c in mine.items() if any(x!=0 for x in c)}
            st = {p:c for p,c in stored_real(k,key).items() if any(x!=0 for x in c)}
            ok = (mine==st) and all(x.is_zero() for x in simple[1]) and all(x.is_zero() for x in simple[2]) and const.is_zero() and all(all(x==0 for x in c[2:]) for c in mine.values())
            tot+=1; bad += (not ok)
            if not ok: print('  MISMATCH k=',k,key)
    print('(key,k) kernels compared',tot,'bad',bad)
