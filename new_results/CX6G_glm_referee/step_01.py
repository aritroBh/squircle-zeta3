#!/usr/bin/env python3
"""Step 01: Verify claim (1) at k=100,200,400"""
import mpmath as mp
mp.mp.dps = 60

def H_mp(s):
    return mp.beta(s/2, mp.mpf(3)/2) / 2

def main_mp_fast(k):
    N = k + 1
    w = mp.e**(mp.pi*1j/4)
    def f(r):
        s = r * w
        return mp.e**(mp.pi*1j*s/2) * (s/(1+s))**N * H_mp(s) / s * w
    return 2/mp.pi * mp.quad(f, [0, mp.inf]).imag

# Load M_k values
table = {}
with open('<project>/research/08_theory/D13f_sweep/mk_deg9.tsv') as f:
    for line in f:
        fields = line.strip().split('\t')
        if fields[0].isdigit():
            import re
            mid_match = re.match(r'\[([-0-9.eE+]+)', fields[1])
            if mid_match:
                table[int(fields[0])] = mp.mpf(mid_match.group(1))

print("Claim (1) verification:")
for k in [100, 200, 400]:
    N = k + 1
    m_k = table[k]
    m_main = main_mp_fast(k)
    rest = abs(m_k - m_main)
    bound = 10600000 * mp.e**(-2 * mp.sqrt(N))
    print(f"k={k}: rest={mp.nstr(rest, 12)}, bound={mp.nstr(bound, 12)}, ratio={mp.nstr(rest/bound, 10)}")
