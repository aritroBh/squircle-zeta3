import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from k_load import *
from math import comb
print("== existence of kernel files k=5..21 and keys")
for k in range(5, 22):
    d = kern(k)
    assert all(key in d for key in KEYS5), k
print("ok")
# (2) no +-i simple pole, zero const
print("== (2) ipos/ineg/const of degree-5 kernels, k=5..21")
for k in range(5, 22):
    row=[]
    for key in KEYS5:
        e = kern(k)[key]
        ip = [x for x in e['ipos']]; ing=[x for x in e['ineg']]
        nz_ip = any(any(c != 0 for c in (x if isinstance(x,(tuple,list)) else (x,))) for x in ip)
        nz_in = any(any(c != 0 for c in (x if isinstance(x,(tuple,list)) else (x,))) for x in ing)
        nz_c = any(c != 0 for c in e['const'])
        row.append((len(ip), len(ing), nz_ip, nz_in, nz_c))
    print(k, row)
