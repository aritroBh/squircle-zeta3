import mpmath as mp, time
mp.mp.dps = 50
t0 = time.time()

def phi(s, c):
    return 1 - (s + c)/mp.sqrt(1 + 2*c*s + s**2)
def dk(k, c):
    # d_k(c) = (-1)^k (1-c^2) P'_{k-1}(c)/k
    f = lambda x: mp.legendre(k-1, x)
    return (-1)**k * (1-c**2) * mp.diff(f, c) / k

# (1) mass: int_0^inf phi = 1-c ; phi(0)=1-c ; decreasing
print("== Lemma B numeric: mass, bounds, tail remainder ==")
for th in [mp.pi/6, mp.pi/3, mp.pi/2, 2*mp.pi/3, 5*mp.pi/6, 0.9*mp.pi]:
    c = mp.cos(th)
    mass = mp.quad(lambda s: phi(s, c), [0, 1, mp.inf])
    # bound constant: phi(s) <= C/s^2 for s>=1, C=(1+|c|)/(2*(1-|c|))
    C = (1+abs(c))/(2*(1-abs(c)))
    viol = mp.findroot(lambda s: phi(s,c)*s**2 - C, (2,5)) if False else None
    maxratio = max(phi(mp.mpf(10)**(mp.mpf(k)/10), c)*(mp.mpf(10)**(mp.mpf(k)/10))**2 for k in range(0,61))
    print(f"theta={mp.nstr(th,6)}: mass-（1-c)={mp.nstr(mass-(1-c),3)}  C_bound={mp.nstr(C,4)}  max s^2 phi(s) on [1,1e6]={mp.nstr(maxratio,4)}  <= C: {maxratio<=C}")

# (2) tail remainder: |phi(s) - sum_{k=2}^L d_k s^{-k}| <= C_L s^{-L-1} for s>=2
print("\ntail remainder constants C_L (sup over s>=2 of s^{L+1}|rem|):")
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
    c = mp.cos(th)
    for L in [2,3,4,6]:
        d = [dk(k,c) for k in range(2, L+1)]
        def rem(s):
            return s**(L+1)*abs(phi(s,c) - sum(d[k-2]*s**(-k) for k in range(2,L+1)))
        vals = [rem(mp.mpf(2)**(mp.mpf(j)/10)) for j in range(0, 121)]
        print(f"  theta={mp.nstr(th,6)} L={L}: sup s>=2 = {mp.nstr(max(vals),4)}")
print("elapsed", time.time()-t0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} Lemma B numeric: mass=1-c (1e-48), explicit s^-2 bound, tail remainder constants  [NUMERICAL]\n")
