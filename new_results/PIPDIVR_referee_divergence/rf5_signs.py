# rf5: sign pattern, odd-k lower bound, and |N_k|^{1/k} vs k/(e pi) (own code).
import sys, re, mpmath as mp
mp.mp.dps = 60
N = {}
for line in open(sys.argv[1]):
    k, r = line.split(" ", 1); N[int(k)] = mp.mpf(re.match(r"\[(\S+) \+/- (\S+)\]", r.strip()).group(1))
pat = {1: 1, 2: 1, 3: -1, 0: -1}
bad = [k for k in range(2, max(N)+1) if mp.sign(N[k]) != pat[k % 4]]
print("k >= 2 violating + + - - (k = 1,2,3,0 mod 4):", bad)
print("odd k with |N_k| < (k-2)!/pi^{k-1}:", [k for k in range(3, max(N)+1, 2) if abs(N[k]) < mp.factorial(k-2)/mp.pi**(k-1)])
for k in (20, 50, 100, 149, 150):
    print("k=%d  |N_k|^{1/k} / (k/(e pi)) = %s" % (k, mp.nstr(abs(N[k])**(mp.mpf(1)/k)/(k/(mp.e*mp.pi)), 6)))
