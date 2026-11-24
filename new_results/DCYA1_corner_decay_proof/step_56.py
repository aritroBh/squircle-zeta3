import sys
sys.path.insert(0,'<project>/research/08_theory/DCYA1_corner_decay_proof/')
import mpmath as mp
mp.mp.dps=50
# Main-term integral on the ray x = e^{i pi/4} t:
# J_j = int e^{-x} x^{N-1} (x+i a)^{-(N-j-3/2)} dx, saddle t_s near sqrt(N a)
# RePhi(t) = -t/sqrt2 + (N-1) log t - (N-j-3/2)/2 * log(t^2 + sqrt2 a t + a^2)
def RePhi(t,N,a,j):
    return -t/mp.sqrt(2)+(N-1)*mp.log(t)-(N-j-mp.mpf(3)/2)/2*mp.log(t**2+mp.sqrt(2)*a*t+a**2)
def dRePhi(t,N,a,j):
    return -1/mp.sqrt(2)+(N-1)/t-(N-j-mp.mpf(3)/2)/2*(2*t+mp.sqrt(2)*a)/(t**2+mp.sqrt(2)*a*t+a**2)
def d2RePhi(t,N,a,j):
    q=t**2+mp.sqrt(2)*a*t+a**2
    return -(N-1)/t**2-(N-j-mp.mpf(3)/2)/2*(2*q-(2*t+mp.sqrt(2)*a)**2)/q**2
for a in [mp.pi/3, 2*mp.pi/3, 3*mp.pi/4, mp.pi/2]:
    N=400; j=0
    # find saddle: Newton on dRePhi
    ts=mp.sqrt(N*a)
    for _ in range(50):
        ts=ts-dRePhi(ts,N,a,j)/d2RePhi(ts,N,a,j)
    curv=-d2RePhi(ts,N,a,j)
    print(f"a={mp.nstr(a,5)}: t_s={mp.nstr(ts,10)} (sqrt(Na)={mp.nstr(mp.sqrt(N*a),8)}), -RePhi''={mp.nstr(curv,8)}, sigma={mp.nstr(1/mp.sqrt(curv),6)}")
    print(f"    RePhi(ts)+theta={mp.nstr(RePhi(ts,N,a,j)+mp.sqrt(2*a*N),8)} (should be O(1) phase decay)")
