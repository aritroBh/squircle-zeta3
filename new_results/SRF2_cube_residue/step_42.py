# Richardson-extrapolate the reliable -qJ values to eta=0: confirm R = I-2 to several digits.
import mpmath as mp
mp.mp.dps = 60
data = [  # (eta, -qJ) from the reliable split-integration
 ('0.2',  '-1.442343324496199313339'),
 ('0.1',  '-1.631773624071964230459'),
 ('0.05', '-1.736024747905954742668'),
 ('0.02', '-1.801875749693849217945'),
 ('0.01', '-1.824398214793221550477'),
 ('0.005','-1.835768857874691694542'),
 ('0.002','-1.842626532892758208994'),
]
pts = [(mp.mpf(a), mp.mpf(b)) for a,b in data]
# The boundary layer makes -qJ(-1+eta) approach R with corrections O(eta^alpha log...).
# Fit R + c1*eta + c2*eta^2 ... on the smallest etas (least boundary-layer curvature is large though).
# Instead fit against eta*(log eta) model? Simplest: polynomial in eta on last 4 points.
import numpy as np
for n, idxs in [(2, [3,4,5]), (3, [3,4,5,6]), (3, [2,3,4,5]), (4, [1,2,3,4,5])]:
    xs = [pts[i][0] for i in idxs]; ys = [pts[i][1] for i in idxs]
    # least-squares polynomial fit degree n-1... use extrapolation at 0 via solve
    A = mp.matrix(len(idxs), n)
    for r, i in enumerate(idxs):
        for c in range(n):
            A[r,c] = pts[i][0]**c
    b = mp.matrix([pts[i][1] for i in idxs])
    coef = mp.lu_solve(A, b)
    print(f"deg-{n-1} fit on etas {[str(pts[i][0]) for i in idxs]}: R = {mp.nstr(coef[0], 18)}")
R = -1-2*mp.pi**mp.mpf('1.5')/mp.gamma(mp.mpf('0.25'))**2
print("exact R =", mp.nstr(R, 20))
