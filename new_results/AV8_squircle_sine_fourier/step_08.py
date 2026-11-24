from mpmath import mp, mpf, pi, sin, cos, quad, matrix, lu_solve
import time, json, os
mp.dps = 45

def pi_p(p):
    p = mpf(p); return 2*pi/(p*sin(pi/p))
def F_p(y, p):
    p=mpf(p); y=mpf(y)
    if y == 0: return mpf(0)
    if y == 1: return pi_p(p)/2
    return y*quad(lambda s: (1-(y*s)**p)**(-1/p), [0,1])
def Phi_p(y,p):
    return pi*F_p(y,p)/pi_p(p)
def b_n_quad(n, p):
    return 4/(n*pi)*quad(lambda y: cos(n*Phi_p(y,p)), [0,1])

ns = [1,3,5,7,9,11,13]
h = mpf('0.02')
offsets = list(range(-4,5))
os.makedirs('scratchpad', exist_ok=True)
ckpt = 'scratchpad/bn_near2.json'
data = {}
if os.path.exists(ckpt):
    data = json.load(open(ckpt))
t00=time.time()
for n in ns:
    for j in offsets:
        key = f"{n},{j}"
        if key in data: continue
        e = j*h
        if j == 0:
            val = mpf(1) if n==1 else mpf(0)
        else:
            val = b_n_quad(n, 2+e)
        data[key] = mp.nstr(val, 42)
        json.dump(data, open(ckpt,'w'))
        print(f"n={n} j={j} done {time.time()-t00:.0f}s", flush=True)
print("ALL DONE", time.time()-t00)
with open('progress.log','a') as f:
    f.write("step_05: computed b_n(2+jh), n odd<=13, j=-4..4, h=0.02, 42 digits\n")
