#!/usr/bin/env python3
"""Step 02: Verify claim (2) at k=100,200,400"""
import mpmath as mp
mp.mp.dps = 60

table = {}
with open('<project>/research/08_theory/D13f_sweep/mk_deg9.tsv') as f:
    for line in f:
        fields = line.strip().split('\t')
        if fields[0].isdigit():
            import re
            mid_match = re.match(r'\[([-0-9.eE+]+)', fields[1])
            if mid_match:
                table[int(fields[0])] = mp.mpf(mid_match.group(1))

def main_theorem_value(N):
    return mp.sqrt(mp.pi)/N * mp.e**(-mp.sqrt(mp.pi*N)) * mp.cos(mp.sqrt(mp.pi*N))

print("Claim (2) verification:")
for k in [100, 200, 400]:
    N = k + 1
    m_k = table[k]
    main_val = main_theorem_value(N)
    lhs = abs(m_k + main_val)
    rhs = 1.21e9 * N**(-1.5) * mp.e**(-mp.sqrt(mp.pi*N))
    print(f"k={k}: lhs={mp.nstr(lhs, 15)}, rhs={mp.nstr(rhs, 15)}, ratio={mp.nstr(lhs/rhs, 10)}")
