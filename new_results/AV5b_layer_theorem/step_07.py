import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)

# Paper numbers to match
C2 = "6*sqrt(2)"
print("C_2 = 6 sqrt(2) =", 6*2**0.5)
print("C_3 = -9 sqrt(3)/8 =", -9*3**0.5/8)

# My tent data derivation for regular 2m-gon: H = sin(pi/2m), alpha = tan(pi/2m)
import math
for m in range(2,7):
    H = math.sin(math.pi/(2*m)); al = math.tan(math.pi/(2*m))
    Cm_paper = 48*m*math.sin(math.pi/(2*m))**4*math.cos(math.pi/(2*m))*(3*math.sin(math.pi/(2*m))**2-1)
    # my formula with overshoot softplus rate r (to be determined): C_m = (3mH/2)(1/r + 3 al/r^2)
    # solve: what (r) with sigma=1 matches?
    import scipy.optimize as so
    f = lambda r: (3*m*H/2)*(1/r + 3*al/r**2) - Cm_paper
    try:
        r = so.brentq(f, 0.01, 100)
        print(f"m={m}: H={H:.6f} alpha={al:.6f} C_paper={Cm_paper:.8f} -> matched r={r:.6f}  r/alpha={r/al:.6f} 1/r={1/r:.6f}")
    except Exception as e:
        # try undershoot sigma1=-1 for the s g1 term: C=(3mH/2)(1/r - 3al/r^2)
        f2 = lambda r: (3*m*H/2)*(1/r - 3*al/r**2) - Cm_paper
        rr = []
        for guess in [0.05,0.1,0.2,0.5,1,2,5,10,50]:
            try:
                r = so.brentq(f2, guess, guess*3+0.01); rr.append(r)
            except Exception: pass
        print(f"m={m}: H={H:.6f} alpha={al:.6f} C_paper={Cm_paper:.8f} overshoot-no-match; undershoot roots r={rr} alpha={al:.6f}")
