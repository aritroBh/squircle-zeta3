import mpmath, json
from mpmath import mp, mpf, sqrt, pi, cos, sin, nstr, sech, matrix, lu_solve
mp.dps = 60
# (a) Minkowski sum square + disk: h = |c|+|s| + r ; P = 8 + 2 pi r exactly
print("(a) Minkowski sum: P = 8 + 2 pi r -- exact by Steiner (no free constants, no zeta)")
# (b) sech-profile family: h_w = h_sq + w * sum_j sech((th - corner_j)/w)
def h_sq(th):
    c, s = cos(th), sin(th)
    return abs(c)+abs(s)
def h_sechfam(th, w):
    corners = [pi/4 + j*pi/2 for j in range(4)]
    return h_sq(th) + w*sum(sech((th-cj)/w) for cj in corners)
def P_sech(w):
    return mpmath.quad(lambda th: h_sechfam(th, w), [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 2*pi])
ws = [mpf('0.02')/mpf('1.5')**j for j in range(6)]
Ps = [P_sech(w) for w in ws]
print("(b) sech profile family:")
pred = lambda w: 8 + 4*pi*w**2
for w,P in zip(ws,Ps):
    print(f"  w={nstr(w,10)} P={nstr(P,40)}  pred 8+4pi w^2={nstr(pred(w),40)}  diff={nstr(P-pred(w),8)}")
