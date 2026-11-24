
import mpmath as mp, json, time
WORK='research/08_theory/K1e_hexagon_closed'
mp.mp.dps=100
direct = json.load(open(f"{WORK}/fresh_direct.json"))
formula = json.load(open(f"{WORK}/fresh_formula.json"))
td = 2*sum(mp.mpf(v) for v in direct.values())
tf = 2*sum(mp.mpf(v) for v in formula.values())
print("total P'''(1): formula =", mp.nstr(tf,50))
print("total P'''(1): direct  =", mp.nstr(td,50))
print("total residual:", mp.nstr(abs(td-tf),3))
# also per-edge from the two files
for k in ['1','2','3']:
    print(f"edge {k}: |formula-direct| = {mp.nstr(abs(mp.mpf(formula[k])-mp.mpf(direct[k])),3)}")
with open(f"{WORK}/progress.log",'a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_25: fresh totals: formula P'''(1)={mp.nstr(tf,30)}, direct={mp.nstr(td,30)}, total residual {mp.nstr(abs(td-tf),3)}.\n")
