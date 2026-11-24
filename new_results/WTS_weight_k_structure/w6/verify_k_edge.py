"""verify step: the 'top-weight k = weight-1 breaks the law' pattern at weights 2 and 3, from OWN definition-based M_k.
 weight 2: law W_2 = (pi^2+6 log^2 2)/24 for k>=2; at k=1, W_2(1) = M_1 (no lower weights).
 weight 3: law W_3(k) = W3' - (k-3) S for k>=3; at k=2, W_3(2) = M_2 - W_2(2) (W_2(2) = law, proved)."""
import sys, os, math, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, W)
import verify_lib as VL
import mpmath as mp
TGT = 30
Y = math.log(10) * (TGT + 12) + 25 + 6
VL.set_dps(TGT + int(Y / 2.302585) + 30)
M1 = VL.direct_M(1, TGT, Y=Y); M2 = VL.direct_M(2, TGT, Y=Y)
VL.set_dps(50)
print("mp.mp.dps =", mp.mp.dps)
lam = mp.log(2); z3 = mp.zeta(3); pi2 = mp.pi ** 2
W2 = (pi2 + 6 * lam ** 2) / 24
S = z3 / 8 + pi2 * lam / 24 + lam ** 3 / 12
W33 = mp.mpf(json.load(open(os.path.join(W, "vals_110", "k3_d2.json")))["value"]) if os.path.exists(os.path.join(W, "vals_110", "k3_d2.json")) else mp.mpf(json.load(open(os.path.join(W, "vals", "k3_d2.json")))["value"])
print("k=1: M_1 = W_2(1) =", mp.nstr(M1, 20), " law value", mp.nstr(W2, 20), " difference", mp.nstr(M1 - W2, 8), "(law FAILS at k=1)")
W32 = M2 - W2
print("k=2: W_3(2) = M_2 - W_2 =", mp.nstr(W32, 20), "; law prediction W3'+S =", mp.nstr(W33 + S, 20), "; difference", mp.nstr(W32 - (W33 + S), 8), "(law FAILS at k=2)")
