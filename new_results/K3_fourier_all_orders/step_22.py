import mpmath as mp
mp.mp.dps = 50
pi = mp.pi; s2 = mp.sqrt(2); z = mp.zeta
exec(open("step_19.py").read().split("Dt = {")[0].replace('mp.dps = 50','mp.dps=50'))  # reuse defs phi, d2phi, d4phi, h, h1, lhs, S0,S2,S4
Dt = {
 (2,0): z(2), (3,0): z(3),
 (4,0): -mp.mpf(33)/16*z(4), (4,1): mp.mpf(7)/2*z(4),
 (5,0): -mp.mpf(3)/4*z(5)-mp.mpf(7)/4*z(2)*z(3), (5,1): 2*z(5)+2*z(2)*z(3),
 (6,0): mp.mpf(4039)/256*z(6)-mp.mpf(7)/8*z(3)**2, (6,1): -mp.mpf(733)/16*z(6)+z(3)**2, (6,2): mp.mpf(93)/2*z(6),
 (7,0): mp.mpf(129)/16*z(7)+mp.mpf(141)/16*z(2)*z(5)+mp.mpf(1071)/160*z(2)**2*z(3),
 (7,1): -23*z(7)-25*z(2)*z(5)-mp.mpf(181)/10*z(2)**2*z(3),
 (7,2): 24*z(7)+24*z(2)*z(5)+mp.mpf(84)/5*z(2)**2*z(3),
}
def comb(tau, order):
    tot = mp.mpf(0)
    for m in range(4):
        for n in range(2*m+2, order+1):
            tot += (pi/4)*(s2/pi)*mp.mpf(2)**(-2*m)/mp.factorial(2*m)*Ss[m]*Dt[(n,m)]*tau**n
    return tot

# predicted limit of [lhs - comb(N)]/tau^(N+1)
def predlim(N):
    n = N+1
    return (pi/4)*(s2/pi)*mp.fsum(mp.mpf(2)**(-2*m)/mp.factorial(2*m)*Ss[m]*Dt[(n,m)] for m in range((n-2)//2+1))

taus = [mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005')]
print("N, residuals at tau=.02,.01,.005, Richardson-extrapolated limit, predicted limit")
for N in [3, 4, 5, 6]:
    rs = [(lhs(t) - comb(t, N))/t**(N+1) for t in taus]
    # linear-in-tau Richardson: r(tau) = L + c tau + ...
    L1 = 2*rs[1]-rs[0]; L2 = 2*rs[2]-rs[1]; L3 = (4*L2-L1)/3
    print(f"N={N}: {[mp.nstr(r,8) for r in rs]}  extrap={mp.nstr(L3,10)}  pred={mp.nstr(predlim(N),10)}  diff={mp.nstr(abs(L3-predlim(N)),3)}")
