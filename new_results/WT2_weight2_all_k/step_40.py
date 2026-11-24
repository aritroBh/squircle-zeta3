src = open('buildintegrand.py').read()
old = """def sqrt_part_PF(k, S):
    \"\"\"PF of  -sqrt2 8^k 4^{-D} v^{2k-1-D} Stilde(v) (v^2+1)^{1-2k} / ((v^2-a^2)(v^2-b^2)).\"\"\"
    St, D = Stilde(S)
    if not St:   # S == 0
        return pf_zero()
    e = 2*k - 1 - D
    num = qsi_pscale(St, SQ2Q * QSI(Fr(-(8**k), 4**D)))
    den_extra_v = 0
    if e >= 0:
        for _ in range(e):
            num = qsi_pmul(num, VPOLY)
    else:
        den_extra_v = -e
    den = qsi_pmul(qsi_pmul(VMINUSA, VPLUSA), qsi_pmul(VMINUSB, VPLUSB))
    for _ in range(den_extra_v):
        den = qsi_pmul(den, VPOLY)
    for _ in range(2*k-1):
        den = qsi_pmul(den, V2P1)
    return pf_from_poly(num, den)"""
new = """def sqrt_part_PF(k, S):
    \"\"\"PF of  -(sqrt2/2) 8^k 4^{-D} v^{2k-1-D} Stilde(v) (v^2+1)^{2-2k} / ((v^2-a^2)(v^2-b^2)).\"\"\"
    St, D = Stilde(S)
    if not St:   # S == 0
        return pf_zero()
    e = 2*k - 1 - D
    num = qsi_pscale(St, SQ2Q * QSI(Fr(-(8**k), 2*4**D)))
    den_extra_v = 0
    if e >= 0:
        for _ in range(e):
            num = qsi_pmul(num, VPOLY)
    else:
        den_extra_v = -e
    den = qsi_pmul(qsi_pmul(VMINUSA, VPLUSA), qsi_pmul(VMINUSB, VPLUSB))
    for _ in range(den_extra_v):
        den = qsi_pmul(den, VPOLY)
    for _ in range(2*k-2):
        den = qsi_pmul(den, V2P1)
    return pf_from_poly(num, den)"""
assert old in src
open('buildintegrand.py','w').write(src.replace(old, new))
print("patched")
