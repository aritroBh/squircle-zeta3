import sys
sys.dont_write_bytecode = True
W = '<project>/research/08_theory/WT3_weight3_all_k'
sys.path.insert(0, '<project>/paper/squircle_note/scripts'); sys.path.insert(0, W)
from fractions import Fraction as Fr
import fastred
from fastred import (QSI, ZERO, AQ, BQ, pf_zero, pf_add, pf_sub, pf_scale, pf_split, pf_mul_DL, pf_eval_at, DL1, DL2, pf_neg)
from buildintegrand import build_integrand_pf
def pf_is_zero(P): return not P[0] and all(e.is_zero() for e in P[1]) and all(e.is_zero() for e in P[2]) and P[3].is_zero()
def simple_part(P): return pf_split(P)[0]
def gms_at(k):
    work = {}
    for key, R in build_integrand_pf(k).items():
        kk, RR = fastred.fold_key(key, R); work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    gms = {}
    for deg in range(k, -1, -1):
        for key in sorted([kk for kk in work if sum(kk)==deg], reverse=True):
            R = work.pop(key)
            if deg == 0: continue
            simple, higher, const = pf_split(R)
            if not higher[1] and not higher[2]: continue
            G = fastred.pf_integrate_higher(higher); C = pf_eval_at(G, BQ); Gm = [G[0], G[1], G[2], -C]; gms[key] = Gm
            i, j = key
            if i >= 1:
                kk, RR = fastred.fold_key((i-1,j), pf_scale(pf_mul_DL(Gm, DL1), Fr(-i))); work[kk] = pf_add(work.get(kk, pf_zero()), RR)
            if j >= 1:
                kk, RR = fastred.fold_key((i,j-1), pf_scale(pf_mul_DL(Gm, DL2), Fr(-j))); work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    return gms
g3, g7 = gms_at(3), gms_at(7)
G3, G21 = g3[(3,0)], g3[(2,1)]
s30 = [{}, [QSI(0,Fr(1,12))], [QSI(0,Fr(1,12))], QSI(Fr(-1,12))]
s21 = [{}, [QSI(0,Fr(-1,4))], [QSI(0,Fr(-1,4))], QSI(Fr(1,4))]
R20_3 = pf_zero(); R20_3[0][ZERO] = QSI(Fr(1,4),Fr(-1,4)); R20_3[0][AQ] = QSI(Fr(1,4)); R20_3[0][-AQ] = QSI(Fr(-1,4)); R20_3[0][-BQ] = QSI(Fr(-1,2))
R11_3 = pf_zero(); R11_3[0][ZERO] = QSI(Fr(1,4),Fr(1,4)); R11_3[0][-AQ] = QSI(Fr(-1,2))
S20 = pf_zero(); S20[0][-AQ] = QSI(Fr(1,2)); S20[0][-BQ] = QSI(Fr(-1,2))
S11 = pf_zero(); S11[0][ZERO] = QSI(Fr(1,2)); S11[0][-AQ] = QSI(-1)
C1_bug = pf_is_zero(pf_sub(simple_part(pf_add(pf_scale(pf_mul_DL(G3, DL1), Fr(-3)), pf_neg(pf_mul_DL(g7[(2,1)], DL2)))), R20_3))
C1 = pf_is_zero(pf_sub(simple_part(pf_add(pf_scale(pf_mul_DL(G3, DL1), Fr(-3)), pf_neg(pf_mul_DL(G21, DL2)))), R20_3))
C2 = pf_is_zero(pf_sub(simple_part(pf_add(pf_scale(pf_mul_DL(s30, DL1), Fr(-3)), pf_neg(pf_mul_DL(s21, DL2)))), S20))
C3 = pf_is_zero(pf_sub(simple_part(pf_scale(pf_mul_DL(G21, DL1), Fr(-2))), R11_3))
C4 = pf_is_zero(pf_sub(simple_part(pf_scale(pf_mul_DL(s21, DL1), Fr(-2))), S11))
print('C1 as in verifier (Gm21 from k=7):', C1_bug)
print('C1 with Gm21^(3) (as in proof_wt3.md):', C1, '| C2:', C2, '| C3:', C3, '| C4:', C4)
# affine Gm law consistency at k=7: Gm^(7) = Gm^(3) + 4 s
aff30 = pf_is_zero(pf_sub(g7[(3,0)], pf_add(G3, pf_scale(s30, Fr(4)))))
aff21 = pf_is_zero(pf_sub(g7[(2,1)], pf_add(G21, pf_scale(s21, Fr(4)))))
print('affine Gm law at k=7 (Gm30, Gm21):', aff30, aff21)
