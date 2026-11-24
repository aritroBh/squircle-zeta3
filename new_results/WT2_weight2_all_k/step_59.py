src = open('fastred.py').read()
old = """    N = ptrimz(list(N)); D = ptrimz(list(D))
    assert len(D) >= 1
    assert len(N) < len(D), ("improper fraction", len(N), len(D))
    P = pf_zero()"""
new = """    N = ptrimz(list(N)); D = ptrimz(list(D))
    assert len(D) >= 1
    P = pf_zero()
    if len(N) >= len(D):
        # polynomial division: keep only constant quotient term
        q = [ZERO]*(len(N)-len(D)+1)
        rem = list(N)
        for t in range(len(N)-len(D), -1, -1):
            c = rem[len(D)-1+t] / D[-1]
            q[t] = c
            for s in range(len(D)):
                rem[s+t] = rem[s+t] - c*D[s]
        while rem and rem[-1].is_zero(): rem.pop()
        assert len(q) <= 1, ("nonconstant polynomial part", q)
        if q: P[3] = q[0]
        N = rem"""
assert old in src
open('fastred.py','w').write(src.replace(old, new))
print("patched pf_from_poly for polynomial part")
