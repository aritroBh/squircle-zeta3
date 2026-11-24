src = open('buildintegrand.py').read()
old = """    for i in range(k, (k-1)//2, -1):
        j = k - i
        P = sqrt_part_PF(k, S[i])            # unfolded l1^i l2^j
        Sj = p_comp_1mt(S[j])                # folded l1^j l2^i  (t -> 1-t)
        P = pf_add(P, sqrt_part_PF(k, Sj))
        if j == 0:                           # key (k,0): A-part + folded B-part
            P = pf_add(P, AB_part_PF(k))
            P = pf_add(P, AB_part_PF(k))
        out[(i, j)] = P"""
new = """    for i in range(k, (k-1)//2, -1):
        j = k - i
        P = sqrt_part_PF(k, S[i])            # unfolded l1^i l2^j
        if i != j:
            Sj = p_comp_1mt(S[j])            # folded l1^j l2^i  (t -> 1-t)
            P = pf_add(P, sqrt_part_PF(k, Sj))
        if j == 0:                           # key (k,0): A-part + folded B-part
            P = pf_add(P, AB_part_PF(k))
            P = pf_add(P, AB_part_PF(k))
        out[(i, j)] = P"""
assert old in src
open('buildintegrand.py','w').write(src.replace(old, new))
print("patched")
