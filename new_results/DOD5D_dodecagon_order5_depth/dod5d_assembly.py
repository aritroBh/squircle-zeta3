"""dod5d_assembly.py -- assemble P^(5)(1)/u two ways and compare to the engine value.
(a) from itab monomial integrals (direct tanh-sinh quadrature);
(b) from the exact G-rewrite + numerical G-table (independent evaluation).
Also assembles the weight-5 PART (for reference) and splits blocks.
"""
import pickle, sys, time
sys.path.insert(0, ".")
import mpmath as mp
import dod5d_grewrite as GR

DPS = 130
mp.mp.dps = DPS + 20
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
rewrite = pickle.load(open("rewrite_exact.pkl","rb"))
r3 = mp.sqrt(3)
u = mp.sqrt(2-r3)
CA = mp.log((2-r3)/2); CC0 = mp.log(2-r3); CC1 = mp.log((2*r3-3)/2)
P5_ENGINE = mp.mpf("-62.53724775498639535309129282205803572684769277410373262688145818074695731422874528214220721559218128482460652943087962941126065293896549327018536685224048730049659698854")
P5u = P5_ENGINE/u

def q3f(c):
    return mp.mpf(str(c[0])) + mp.mpf(str(c[1]))*r3

SLOT = {GR.L0: 0, GR.L1: 1, GR.P0: 2, GR.P1: 3, GR.P2: 4, GR.P3: 5}

mode = sys.argv[1]
if mode == 'itab':
    IT = pickle.load(open("itab140.pkl","rb"))
    IT = {m: mp.mpf(v) for m, v in IT.items()}
    ival = IT
else:
    G1raw = pickle.load(open("gtab6_d130_N320.pkl","rb"))
    G1 = {k: mp.mpf(v) for k, v in G1raw.items()}
    ival = {}
    for m, rw in rewrite.items():
        tot = mp.mpf(0)
        for w, cw in rw.items():
            tot += q3f(cw)*G1[tuple(SLOT[x] for x in w)]
        ival[m] = tot

# assemble full P5/u and weight-5 part
full = mp.mpf(0); w5 = mp.mpf(0)
blocks = {}
for (n, le), coef in E5X.items():
    v = q3f(coef)*ival[(n, le[0], le[1], le[2], le[3], le[4], le[5])]
    if le[6]: v *= CA**le[6]
    if le[7]: v *= CC0**le[7]
    if le[8]: v *= CC1**le[8]
    full += v
    if sum(le) == 5:
        w5 += v
        nz = tuple(i for i in range(4) if le[2+i] > 0)
        blk = 'local' if not nz else ('single' if len(nz)==1 else ('double' if len(nz)==2 else ('triple' if len(nz)==3 else 'quad')))
        blocks[blk] = blocks.get(blk, mp.mpf(0)) + v
res = 1440*full - P5u
print("mode:", mode)
print("P5/u assembled:", mp.nstr(1440*full, 60))
print("P5/u engine    :", mp.nstr(P5u, 60))
print("residual:", mp.nstr(abs(res), 3))
print("weight-5 part of P5/u:", mp.nstr(1440*w5, 60))
for k in sorted(blocks):
    print("  block %-7s %s" % (k, mp.nstr(1440*blocks[k], 50)))
