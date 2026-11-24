code = r'''import mpmath as mp, time
import cx6r_common as C
mp.mp.dps = 100
print("CX6R step_10: bound (12) |A_k(s)| <= 2 rho^k/|s+2|^{k+1} + 2 sum_{i>=1} (2^i/i!) rho^{k-i}/|s+i+1|^{k-i+1}")
def RHS12(k, s):
    rho = abs(s)
    out = 2*rho**k/abs(s+2)**(k+1)
    for i in range(1, k+1):
        out += 2*(mp.mpf(2)**i/mp.factorial(i))*rho**(k-i)/abs(s+i+1)**(k-i+1)
    return out
worst = (mp.mpf(0), None)
t=time.time()
for k in [0,1,2,3,5,8,12,16,20,24,28,30]:
    for r in [0, mp.mpf('0.25'), mp.mpf('0.5'), 1, 2, 4, 7, 12, 20, 35, 60]:
        s = mp.mpf('0.5')+r*C.w4
        av = abs(C.Ak(k, s, level=6))
        rhs = RHS12(k, s)
        ratio = av/rhs
        if ratio > worst[0]: worst = (ratio, (k, r, av, rhs))
print(f"({time.time()-t:.0f}s) worst ratio |A_k|/RHS(12) = {mp.nstr(worst[0],10)} at k={worst[1][0]}, r={mp.nstr(worst[1][1],6)}")
print(f"   |A_k| = {mp.nstr(worst[1][2],10)}, RHS = {mp.nstr(worst[1][3],10)}")
# second-worst region scan near the worst
k0, r0 = worst[1][0], worst[1][1]
w2 = mp.mpf(0)
for k in range(max(0,k0-3), k0+4):
    for j in range(41):
        r = r0*(mp.mpf('0.5') + j*mp.mpf('0.025'))
        s = mp.mpf('0.5')+r*C.w4
        ratio = abs(C.Ak(k, s, level=6))/RHS12(k, s)
        w2 = max(w2, ratio)
print("refined worst ratio near worst point:", mp.nstr(w2,10))
'''
open('step_10.py','w').write(code)
exec(compile(code,'step_10.py','exec'))
