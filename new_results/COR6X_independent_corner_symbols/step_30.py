t = r'''
import mpmath as mp
mp.mp.dps = 50
def I2(a0, a1, a2, a3):
    return mp.quad(lambda t: mp.log((t - a2)/(a0 - a2)) / (t - a1), [a0, a3])
def I1(a0, a1, a2):
    return mp.log((a2 - a1)/(a0 - a1))
a0, a1, a2, a3 = mp.mpf(0), mp.mpf('-0.31'), mp.mpf('-0.17'), mp.mpf('0.8')
val = I2(a0, a1, a2, a3)
for k in range(4):
    pp = [a0, a1, a2, a3]
    print(f"--- parameter a{k}")
    for he in ('1e-6','1e-8','1e-10'):
        h = mp.mpf(he)
        pp2 = pp[:]; pp2[k] += h
        dnum = (I2(*pp2) - val) / h
        print(f"   h={he}: {mp.nstr(dnum, 15)}")
# exact formula predictions
prev = {1: 3, 2: 1}; nxt = {1: 2, 2: 0}
omit = {1: I1(a0, a2, a3), 2: I1(a0, a1, a3)}
params = [a0, a1, a2, a3]
for k in range(4):
    pred = mp.mpf(0)
    for i in (1, 2):
        term = mp.mpf(0)
        if k == prev[i]: term += 1/(params[prev[i]] - params[i])
        if k == i:       term += -1/(params[prev[i]] - params[i]) + 1/(params[i] - params[nxt[i]])
        if k == nxt[i]:  term += -1/(params[i] - params[nxt[i]])
        pred += omit[i] * term
    print(f"formula d/da{k} = {mp.nstr(pred, 15)}")
'''
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","-c",t],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=300)
print(r.stdout); print("ERR:", (r.stderr or "")[-800:])
