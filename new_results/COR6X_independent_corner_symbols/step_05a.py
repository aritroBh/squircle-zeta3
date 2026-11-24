
# step_05a.py -- numerical verification of the differential formula (weight 2 + one weight-3)
import mpmath as mp
mp.mp.dps = 40

def I2(a0, a1, a2, a3):
    return mp.quad(lambda t: mp.log((t - a2)/(a0 - a2)) / (t - a1), [a0, a3])
def I1(a0, a1, a2):
    return mp.log((a2 - a1)/(a0 - a1))

# all letters off the path [0, 0.8]
a0, a1, a2, a3 = mp.mpf(0), mp.mpf('-0.31'), mp.mpf('-0.17'), mp.mpf('0.8')
params = [a0, a1, a2, a3]
val = I2(*params)

h = mp.mpf('1e-16')
# cyclic recursion data: for i=1..2: (prev_i, next_i); prev_1 = a3, next_2 = a0
prev = {1: 3, 2: 1}
nxt = {1: 2, 2: 0}
omit = {1: lambda: I1(a0, a2, a3), 2: lambda: I1(a0, a1, a3)}
for k in range(4):
    pp = params[:]; pp[k] += h
    dnum = (I2(*pp) - val) / h
    pred = mp.mpf(0)
    for i in (1, 2):
        # d/dp_k [ log(prev_i - a_i) - log(next_i - a_i) ]
        term = mp.mpf(0)
        if k == prev[i]: term += 1/(params[prev[i]] - params[i])
        if k == i:       term += -1/(params[prev[i]] - params[i]) - 1/(params[i] - params[nxt[i]])
        if k == nxt[i]:  term += +1/(params[i] - params[nxt[i]])
        pred += omit[i]() * term
    print(f"d/da{k}: numeric {mp.nstr(dnum,10)}  formula {mp.nstr(pred,10)}  err {mp.nstr(abs(dnum-pred),3)}")

# weight 3 spot check: vary only a1 and a4 (endpoint)
def I3_inner(a0, a2, a3, t1):
    return mp.quad(lambda t2: mp.log((t2 - a3)/(a0 - a3)) / (t2 - a2), [a0, t1])
def I3(a0, a1, a2, a3, a4):
    return mp.quad(lambda t1: I3_inner(a0, a2, a3, t1) / (t1 - a1), [a0, a4])
b0, b1, b2, b3, b4 = mp.mpf(0), mp.mpf('-0.31'), mp.mpf('-0.17'), mp.mpf('-0.05'), mp.mpf('0.8')
print("weight-3 check:")
for k, name, hfun in ((1, 'a1', lambda: I3(b0, b1 + h, b2, b3, b4)),
                      (4, 'a4', lambda: I3(b0, b1, b2, b3, b4 + h)),
                      (0, 'a0', lambda: I3(b0 + h, b1, b2, b3, b4)),
                      (2, 'a2', lambda: I3(b0, b1, b2 + h, b3, b4)),
                      (3, 'a3', lambda: I3(b0, b1, b2, b3 + h, b4))):
    dnum = (hfun() - I3(b0, b1, b2, b3, b4)) / h
    # cyclic: prev_1 = a4; next_3 = a0
    prevv = {1: 4, 2: 1, 3: 2}; nxtt = {1: 2, 2: 3, 3: 0}
    pp = [b0, b1, b2, b3, b4]
    pred = mp.mpf(0)
    for i in (1, 2, 3):
        Iom = I2(b0, *[pp[j] for j in (1,2,3) if j != i], b4)
        term = mp.mpf(0)
        if k == prevv[i]: term += 1/(pp[prevv[i]] - pp[i])
        if k == i:        term += -1/(pp[prevv[i]] - pp[i]) - 1/(pp[i] - pp[nxtt[i]])
        if k == nxtt[i]:  term += +1/(pp[i] - pp[nxtt[i]])
        pred += Iom * term
    print(f"  d/d{name}: numeric {mp.nstr(dnum,10)}  formula {mp.nstr(pred,10)}  err {mp.nstr(abs(dnum-pred),3)}")
