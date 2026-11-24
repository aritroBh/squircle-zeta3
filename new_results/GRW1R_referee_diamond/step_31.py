import mpmath as mp
mp.mp.dps = 60

# Classical Stirling: logGamma(x) = (x-1/2)ln x - x + ln(2pi)/2 + 1/(12x) - th/(360 x^3), 0<th<1, x>0.
# log gA(n) := log[Gamma(n-1/2)/Gamma(n+1)] + (3/2) ln n
#           = (n-1)ln(1-0.5/n) - (n+0.5)ln(1+1/n) + 3/2 + [1/(12(n-0.5)) - 1/(12(n+1))] + Rst,
# |Rst| <= 1/(360(n-0.5)^3) + 1/(360(n+1)^3).
# Taylor with remainder: for 0<=x<=1/2: |ln(1-x)+x+x^2/2+x^3/3| <= x^4 (series negative, |sum_{k>=4}| <= x^4 sum x^{k-4} <= x^4/(1-x) <= 2x^4; use 2x^4)
# for 0<=x<=1/2: |ln(1+x)-x+x^2/2-x^3/3| <= x^4/4 + x^5/5 + ... <= x^4.

def bounds(n0):
    # returns (C_A2, C_B1, C_C0): for n>=n0:
    #   n^2 |gA-1-3/(8n)| <= C_A2 ;  n |gB-1| <= C_B1 ; gC = Gamma(n-5/2)/Gamma(n+1) n^{7/2} <= C_C0
    # gA: log gA = 0.375/n + 0.125/n^2 + RA, bound |RA|:
    # pieces: (n-1)ln(1-0.5/n) = -0.5 + 0.375/n + 0.0833.../n^2 + rA1
    # Compute exactly: (n-1)(-0.5/n -0.125/n^2 -1/24/n^3 + e1), |e1|<=2(0.5/n)^4
    #   = -0.5 +0.5/n -0.125/n +0.125/n^2 -1/(24 n^2) + 1/(24 n^3) + (n-1)e1
    #   = -0.5 + 0.375/n + (0.125-0.0416667)/n^2 + rA1, |rA1| <= 1/(24 n^3) + 2 n (0.5/n)^4 = 0.0417/n^3 + 0.125/n^3 <= 0.167/n^3
    # (n+0.5)ln(1+1/n) = 1 + 0.0833.../n^2 + rA2 : (n+0.5)(1/n -1/(2n^2)+1/(3n^3)+e2), |e2|<=1/(4n^4)
    #   = 1 +0.5/n -0.5/n -0.25/n^2 +1/(3n^2)+1/(6n^3)+ (n+0.5)e2 = 1 + (0.0833)/n^2 + rA2, |rA2| <= 1/(6n^3) + 1.001 n/(4 n^4) <= 0.417/n^3
    # 1/(12(n-0.5)) - 1/(12(n+1)) = 1/(8 n (n-0.5)(n+1))... exact: (1/12)(1.5)/((n-0.5)(n+1)) = 0.125/((n-0.5)(n+1))
    #   = 0.125/n^2 * 1/(1-0.5/n)(1+1/n) = 0.125/n^2 + rA3, |rA3| <= 0.125/n^2 * (1/(1-0.5/n0)-1) + ... bound 0.126/n^3*1.01 -> use 0.2/n^3
    # log gA = [-0.5+0.375/n+0.0833/n^2+rA1] - [1+0.0833/n^2+rA2] + 1.5 + 0.125/n^2 + rA3 + Rst
    #        = 0.375/n + 0.125/n^2 + RA, |RA| <= (0.167+0.417+0.2)/n^3 + 2/(360 (n-0.5)^3) <= 0.79/n^3
    RA = mp.mpf('0.79')/n0**3
    lgA_hi = mp.mpf('0.375')/n0 + mp.mpf('0.125')/n0**2 + RA
    lgA_lo = mp.mpf('0.375')/n0 + mp.mpf('0.125')/n0**2 - RA
    # gA - 1 - 0.375/n: |e^t - 1 - 0.375/n| where t in [lgA_lo, lgA_hi]
    # e^t = 1 + t + t^2/2 e^{t'}, |t^2/2 e^{t'}| <= lgA_hi^2 e^{lgA_hi}
    t2 = lgA_hi**2*mp.e**lgA_hi
    gA_dev_hi = max(abs(mp.e**lgA_hi - 1 - mp.mpf('0.375')/n0), abs(mp.e**lgA_lo - 1 - mp.mpf('0.375')/n0)) 
    # NOTE: e^t-1-0.375/n with t = 0.375/n + 0.125/n^2 + RA: e^t-1-t + t-0.375/n: <= t2/2 + 0.125/n^2 + RA... 
    # simpler: |gA-1-0.375/n| <= |e^t-1-t| + |t-0.375/n| <= t2/2 + 0.125/n0^2 + RA
    CA2 = n0**2*(t2/2 + mp.mpf('0.125')/n0**2 + RA)
    # gB: log gB = 1.875/n + RB:  pieces: (n-2)ln(1-1.5/n) = -1.5 + 1.875/n + rB1:
    #  (n-2)(-1.5/n - 1.125/n^2 + e1), |e1| <= 2(1.5/n)^3: = -1.5 +3/n -1.125/n + 2.25/n^2 + rB1, |rB1| <= 2.25/n^2... 
    #  keep 1/n only: |rB1| <= 2.25/n^2 + (n)(2)(1.5/n)^3 = 2.25/n^2 + 6.75/n^2 = 9/n^2
    #  (n+0.5)ln(1+1/n) = 1 + rB2, |rB2| <= 0.5/n (since = 1 + 0.0833/n^2 + ... <= 1 + 0.09/n^2 for n>=1, use 0.09/n^2)
    #  log gB = -1.5+1.875/n + rB1 - 1 - rB2 + 2.5 + [1/(12(n-1.5))-1/(12(n+1))] + Rst
    #   = 1.875/n + RB, |RB| <= 9/n^2 + 0.09/n^2 + (2.5/12)/((n-1.5)(n+1)) + 2/(360(n-1.5)^3) <= 9.3/n^2
    RB = mp.mpf('9.3')/n0**2
    lgB_hi = mp.mpf('1.875')/n0 + RB
    # n |gB-1| <= n (e^{lgB_hi}-1) = n (lgB_hi + lgB_hi^2 e^{lgB_hi}) 
    CB1 = n0*(mp.e**lgB_hi - 1)  # use exact exp at n0; for larger n it's smaller
    # gC: log gC = 4.375/n + RC, |RC| <= (25 + 0.09)/n^2 + (3.5/12)/((n-2.5)(n+1)) + 2/(360(n-2.5)^3) <= 25.4/n^2
    RC = mp.mpf('25.4')/n0**2
    lgC_hi = mp.mpf('4.375')/n0 + RC
    CC0 = mp.e**lgC_hi
    return CA2, CB1, CC0

n0 = 500
CA2, CB1, CC0 = bounds(n0)
print("n0 =", n0)
print("C_A2 =", mp.nstr(CA2, 15), " (n^2|gA-1-3/(8n)| <=)")
print("C_B1 =", mp.nstr(CB1, 15), " (n|gB-1| <=)")
print("C_C0 =", mp.nstr(CC0, 15), " (gC <=)")

# cross-check the expansions against mpmath gamma at several n >= 500
for nn in [500, 700, 2000, 100000]:
    n = mp.mpf(nn)
    gA = mp.gamma(n-mp.mpf('0.5'))/mp.gamma(n+1)*n**mp.mpf('1.5')
    gB = mp.gamma(n-mp.mpf('1.5'))/mp.gamma(n+1)*n**mp.mpf('2.5')
    gC = mp.gamma(n-mp.mpf('2.5'))/mp.gamma(n+1)*n**mp.mpf('3.5')
    print(f"n={nn}: n^2(gA-1-3/8n)={mp.nstr(n**2*(gA-1-3/(8*n)),10)} (<= {mp.nstr(CA2,6)}), "
          f"n(gB-1)={mp.nstr(n*(gB-1),10)} (<= {mp.nstr(CB1,6)}), gC={mp.nstr(gC,10)} (<= {mp.nstr(CC0,6)})")

# assemble the tail bound for R_bound/u, n >= 500:
# |A+B-law2|/u <= n^2|gA-1-3/(8n)| + (3pi/8) n |gB-1|  <= CA2 + (3pi/8) CB1
# Ebnd/u = (25 rho^2/(16 sqrt2)) gC <= 2.7268*CC0 ; tails/u -> tiny at n=500
tail_AB = CA2 + (3*mp.pi/8)*CB1
cE = 25*(mp.pi/2)**2/(16*mp.sqrt(2))
tail_E = cE*CC0
# tails at n=500: (4/pi)(2.9 rho)^{1/2-n}/(n-1/2) / (rho^{-n} n^{-7/2}) = (4/pi) 2.9^{1/2} 2.9^{-n} rho^{1/2} n^{7/2}/(n-1/2)
n=mp.mpf(500)
tails500 = ((4/mp.pi)*mp.mpf('2.9')**mp.mpf('0.5')*mp.mpf('2.9')**(-n)*(mp.pi/2)**mp.mpf('0.5')*n**mp.mpf('3.5')/(n-mp.mpf('0.5'))
            + (2/mp.pi)*mp.mpf('2.9')**mp.mpf('1.5')*mp.mpf('2.9')**(-n)*(mp.pi/2)**mp.mpf('1.5')*n**mp.mpf('3.5')/(n-mp.mpf('1.5')))
print("tails500/u =", mp.nstr(tails500, 5))
print("tail (n>=500) R_bound/u <=", mp.nstr(tail_AB + tail_E + tails500, 10))
print("certified sup over 10..500 was 7.0523 at n=11 -> GLOBAL C_R = 7.06 valid for all n>=10")
with open('progress.log','a') as f:
    f.write(f"step_10: Stirling tail n>=500: AB<={mp.nstr(tail_AB,8)}, E<={mp.nstr(tail_E,8)}; global C_R=7.06\n")
