# Referee item 2: Lemma Z on a w-grid. E_w(zeta) = 1 + e^{w zeta} + e^{zeta}.
# For each w: all zeros in |zeta| < RAD by Newton multistart; count checked by the argument principle.
import numpy as np
RHO = np.sqrt(np.pi**2 + np.log(2)**2)
RAD = 7.0

def E(z, w):  return 1 + np.exp(w*z) + np.exp(z)
def dE(z, w): return w*np.exp(w*z) + np.exp(z)

def winding(w, rad=RAD, n=20000):
    th = np.linspace(0, 2*np.pi, n+1)
    v = E(rad*np.exp(1j*th), w)
    return int(round(np.sum(np.diff(np.unwrap(np.angle(v))))/(2*np.pi)))

xs = np.linspace(-RAD, RAD, 41)
STARTS = (xs[:, None] + 1j*xs[None, :]).ravel()
STARTS = STARTS[np.abs(STARTS) < RAD + 1]

def zeros(w):
    z = STARTS.copy()
    for _ in range(80):
        with np.errstate(all='ignore'):
            z = z - E(z, w)/dE(z, w)
    ok = np.isfinite(z) & (np.abs(E(z, w)) < 1e-10) & (np.abs(z) < RAD)
    out = []
    for c in z[ok]:
        if all(abs(c-d) > 1e-7 for d in out): out.append(c)
    return np.array(out)

ws = np.concatenate([np.linspace(0, 1, 2001), [1e-6, 1e-5, 1e-4, 1-1e-4, 1-1e-5, 1-1e-6]])
ws = np.unique(ws)
mismatch = 0; rows = []
for w in ws:
    Z = zeros(w)
    nw = winding(w)
    if nw != len(Z): mismatch += 1
    rows.append((w, np.min(np.abs(Z))/RHO, np.min(np.abs(Z.imag))/np.pi, len(Z), nw))
rows = np.array(rows)
print("w-points:", len(ws), " Newton count != winding count at", mismatch, "points")
print("global min |zeta|/rho* =", rows[:, 1].min(), " at w =", rows[np.argmin(rows[:, 1]), 0])
inner = (rows[:, 0] > 0) & (rows[:, 0] < 1)
print("min over 0<w<1 of |zeta|/rho* - 1 =", (rows[inner, 1]-1).min())
print("min |Im zeta|/pi over all w =", rows[:, 2].min())
for w0 in [0, 1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.19, 0.5, 0.81, 0.95, 0.999, 1]:
    i = np.argmin(abs(rows[:, 0]-w0))
    print(f"  w={rows[i,0]:.6g}  min|zeta|/rho*={rows[i,1]:.9f}  (1+w/2 ~ {1+min(rows[i,0],1-rows[i,0])/2:.6f})  #zeros={int(rows[i,3])}")
# symmetry w -> 1-w
a = rows[(rows[:, 0] <= 0.5) & np.isclose(rows[:, 0]*2000, np.round(rows[:, 0]*2000))]
b = {round(r[0], 6): r[1] for r in rows}
print("max |m(w)-m(1-w)| on grid:", max(abs(b[round(r[0], 6)] - b[round(1-r[0], 6)]) for r in a))
# the case-(ii) inequality a >= log(2 cos beta) at zeros with pi <= |b| < rho*, a >= 0
viol = 0; cnt = 0
for w in ws[::5]:
    for z in zeros(w):
        for zz in (z, -z):  # include the reflected zero used for a<0
            a_, b_ = zz.real, abs(zz.imag)
            if np.pi <= b_ < RHO and a_ >= 0:
                cnt += 1
                if a_ < np.log(2*np.cos(b_-np.pi)) - 1e-12: viol += 1
print("band zeros checked:", cnt, " violations of a>=log(2cos beta):", viol)
