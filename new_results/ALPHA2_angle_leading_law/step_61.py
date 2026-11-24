from mpmath import mp, mpf, nstr, pi, sin, cos, legenp
mp.dps = 55
def Jfull_series3(s, ca, N=3000):
    x = -ca
    P0, P1 = mpf(1), x
    tot = P0*(1/s + 1/(1-s)) + P1*(1/(1+s) + 1/(2-s))
    for n in range(1, N-1):
        P0, P1 = P1, ((2*n+1)*x*P1 - n*P0)/(n+1)
        tot += P1*(1/(n+1+s) + 1/(n+2-s))
    # tail estimate: P_n ~ O(n^-1/2), terms ~ P_n*(2n+3)/((n+s)(n+1-s)) ~ n^-2.5 -> tail ~ N^-1.5
    return tot
for alpha in [pi/2, pi/3, 2*pi/3, pi/6]:
    ca = cos(alpha)
    for s in [mpf('0.9'), mpf('0.3'), mpf('0.5')]:
        a = Jfull_series3(s, ca)
        b = pi*legenp(s-1, 0, -ca, type=2)/sin(pi*s)
        print(f"alpha={nstr(alpha,6)}, s={nstr(s,4)}: |series-cand|={nstr(abs(a-b),3)}")
import subprocess
print(subprocess.run(['tail','-3','<project>/research/08_theory/ALPHA2_angle_leading_law/kern18_old.out'], capture_output=True, text=True).stdout)
