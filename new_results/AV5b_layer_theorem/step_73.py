import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
step = r'''
import mpmath as mp
mp.mp.dps = 55
import json, sys

VS = [(1,1),(1,-1),(-1,1),(-1,-1)]
def mfun(th):
    c, s = mp.cos(th), mp.sin(th)
    return abs(c)+abs(s)
def h1(th, tau):
    c, s = mp.cos(th), mp.sin(th)
    m = abs(c)+abs(s)
    if m == 0: return m
    L = mp.fsum(mp.exp((v[0]*c+v[1]*s-m)/tau) for v in VS)
    return m*mp.exp(tau*mp.log(L))
def h2(th, tau):
    c, s = mp.cos(th), mp.sin(th)
    return tau*mp.log(mp.fsum(mp.exp((v[0]*c+v[1]*s)/tau) for v in VS))

def convexity_check(hfun, tau):
    # r = h + h'' >= 0 on grid
    worst = mp.mpf('1e9'); argw = 0
    N = 400
    for k in range(N):
        th = 2*mp.pi*k/N
        r = hfun(th, tau) + mp.diff(lambda x: hfun(x, tau), th, 2)
        if r < worst: worst, argw = r, th
    return worst, argw

def P_of(hfun, tau):
    # integrate over [0, 2pi), splitting near the 4 corners (axis directions)
    corners = [mp.mpf(0), mp.pi/2, mp.pi, 3*mp.pi/2]
    W = 20*tau
    pts = {mp.mpf(0), 2*mp.pi}
    for cth in corners:
        pts.add((cth-W) % (2*mp.pi)); pts.add((cth+W) % (2*mp.pi))
    pts = sorted(pts)
    tot = mp.mpf(0)
    for a, b in zip(pts[:-1], pts[1:]):
        if b - a < 1e-30: continue
        tot += mp.quad(lambda th: hfun(th, tau), [a, b])
    return tot

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "convex":
        for tstr in ["0.05","0.02","0.01","0.005"]:
            t = mp.mpf(tstr)
            w, a = convexity_check(h1, t)
            print(f"F1 tau={tstr}: min h+h'' = {mp.nstr(w,12)} at theta={mp.nstr(a,8)}")
    else:
        out = []
        for tstr in sys.argv[2].split(","):
            t = mp.mpf(tstr)
            v1 = P_of(h1, t); v2 = P_of(h2, t)
            out.append((tstr, v1, v2))
            print(f"tau={tstr}: P_F1={mp.nstr(v1,35)}  P_F2={mp.nstr(v2,35)}", flush=True)
        json.dump([[t, str(a), str(b)] for t,a,b in out], open("ckpt_newfam.json","w"))
'''
with open("av5b_newfam.py","w") as f: f.write(step)
import subprocess, time
PY = "<project>/.venv/bin/python"
t0=time.time()
r = subprocess.run([PY,"av5b_newfam.py","convex"],capture_output=True,text=True,timeout=600)
print(r.stdout, r.stderr[-300:])
print("elapsed",time.time()-t0)
