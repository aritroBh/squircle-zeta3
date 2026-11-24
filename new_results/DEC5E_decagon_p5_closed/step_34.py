import pickle
import mpmath as mp
mp.mp.dps = 600
W5 = pickle.load(open('wparts_dec_d520.pkl','rb'))
W8 = pickle.load(open('wparts_dec_d800.pkl','rb'))
for k in range(6):
    a = mp.mpf(W5[f'Wt{k}']); b = mp.mpf(W8[f'Wt{k}'])
    print(f"Wt{k}: |520-800| =", mp.nstr(abs(a-b), 3))
    for c in ['local','single-1','single-4','mixed']:
        d = abs(mp.mpf(W5[f'blk{k}'][c]) - mp.mpf(W8[f'blk{k}'][c]))
        if d > mp.mpf(10)**(-450): print("   block", c, mp.nstr(d,3))
open('progress.log','a').write("step_30/31: G-tables d520/N1750 (23.9s, resid 4.1e-522) and d800/N2700 (47.0s, resid 4.6e-802) built+validated; Wt0..Wt5 assembled at both; sum matches P5/L gold 1.1e-163 (gold limit); Wt0 ~ 1e-512/1e-793 (exact-zero candidate); 520-vs-800 cross-agreement below\n")
