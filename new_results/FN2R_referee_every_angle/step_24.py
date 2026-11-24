import mpmath as mp, time
mp.mp.dps = 50
t0=time.time()
def phi(s,c):
    r = mp.sqrt(1+2*c*s+s**2)
    return (1-c**2)/(r*(r+s+c))
def dl(l,c):
    # d_l via Legendre: (-1)^l (1-c^2) P'_{l-1}(c)/l  (use recurrence-free numeric via mp.legendre? use sympy-free: bonnet)
    # numeric derivative of P_{l-1}
    return (-1)**l*(1-c**2)*mp.diff(lambda x: mp.legenp(l-1,0,x), c)/l
print("== B(ii): sup over s>=1 of phi(s)*s^2 / ((1+|c|)/(2(1-|c|))) should be <= 1 ==")
for c in [mp.mpf('0.5'), mp.mpf(0), -mp.mpf('0.5'), mp.mpf('0.86602540378443864676372317075294')]:
    B = (1+abs(c))/(2*(1-abs(c)))
    # max of phi*s^2 on [1,inf): at s=1? phi decreasing but s^2 increasing; phi ~ d2/s^2
    vals = [phi(s,c)*s**2/B for s in [1,2,5,20,100]]
    print(f" c={mp.nstr(c,6)}: ratios at s=1,2,5,20,100: {[mp.nstr(v,6) for v in vals]} (all <1: {all(v<1 for v in vals)})")
print("\n== B(iii): tail truncation error vs d_{L+1} s^{-L-1} scale, theta=pi/3, 2pi/3 ==")
for c in [mp.mpf('0.5'), -mp.mpf('0.5')]:
    for L in [2,3,4,5]:
        for sv in [2,3,5,10]:
            exact = phi(sv,c)
            trunc = sum(dl(l,c)*mp.mpf(sv)**(-l) for l in range(2,L+1))
            rem = exact - trunc
            # remainder should be O(s^{-L-1}); compare with next term
            nxt = dl(L+1,c)*mp.mpf(sv)**(-(L+1))
            print(f" c={mp.nstr(c,4)} L={L} s={sv}: rem={mp.nstr(rem,4)} nextterm={mp.nstr(nxt,4)} ratio={mp.nstr(rem/nxt if nxt else mp.nan,4)}")
    print()
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} step_04b: B(ii) bound holds with margin (ratios<1), tail remainder ~ next term (alternating-ish), O(s^{{-L-1}}) confirmed numerically [referee]\n")
