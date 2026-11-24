wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/verify_full.py') as f:
    content = f.read()

old = """    # -- 0.1 partial_y and partial_x as EXACT finite-sum identities (truncation N).
    N = 25
    S = sum(x**n * y**m / (n**2 * m) for n in range(1, N) for m in range(1, n))
    lhs_y = sp.diff(S, y)
    rhs_y = sp.expand((sum(x**n/n**2 for n in range(1, N))
                       - sum((x*y)**n/n**2 for n in range(1, N))/y)/(1-y))
    r = sp.expand(lhs_y - rhs_y)
    print("  0.1a d_y Li_21 truncated-series identity: residual =", r)
    assert r == 0
    lhs_x = sp.diff(S, x)
    rhs_x = sp.expand(sum(x**n * y**m / (n * m) for n in range(1, N) for m in range(1, n))/x)
    r = sp.expand(lhs_x - rhs_x)
    print("  0.1b d_x Li_21 truncated-series identity: residual =", r)
    assert r == 0
    # -- 0.2 Li_{1,1} closed form:  F = log(1-x)log(1-y) + Li2(-y/(1-y)) - Li2(-y(1-x)/(1-y))
    #     d/dx F = -log(1-x y)/(1-x) (the series derivative of Li_{1,1}); F(0,y)=0.
    z = -y*(1-x)/(1-y)
    dz = sp.diff(z, x)
    dF = (-sp.log(1-y)/(1-x)
          - (-sp.log(1-z)) * dz / z)          # d/dx[-Li2(z)] = +log(1-z) z'/z
    r = sp.factor(dF - (-sp.log(1-x*y)/(1-x)))
    print("  0.2  d/dx Li_11 closed form minus series derivative =", r)
    assert r == 0
    print("       F(0,y) = log(1)*log(1-y) + Li2(-y/(1-y)) - Li2(-y/(1-y)) = 0 (exact)")
    print("PART 0 OK")"""

new = """    # -- 0.1 partial_y and partial_x as EXACT finite-sum identities (truncation N).
    #     d/dy S_N:  (1-y)*d/dy S_N = sum_{n<N} x^n/n^2 (1 - y^{n-1})  (exact poly id);
    #     dividing by (1-y) (formal, y != 1) gives d_y Li_21 = (Li2(x) - Li2(x y)/y)/(1-y).
    N = 25
    S = sum(x**n * y**m / (n**2 * m) for n in range(1, N) for m in range(1, n))
    lhs_y = sp.diff(S, y)
    num_y = (sum(x**n/n**2 for n in range(1, N))
             - sum((x*y)**n/n**2 for n in range(1, N))/y)
    r = sp.expand((1-y)*lhs_y - num_y)
    print("  0.1a (1-y) d_y Li_21 truncated-series identity: residual =", r)
    assert r == 0
    lhs_x = sp.diff(S, x)
    rhs_x = sp.expand(sum(x**n * y**m / (n * m) for n in range(1, N) for m in range(1, n))/x)
    r = sp.expand(lhs_x - rhs_x)
    print("  0.1b d_x Li_21 truncated-series identity: residual =", r)
    assert r == 0
    # -- 0.2 Li_{1,1} closed form:  F = log(1-x)log(1-y) + Li2(-y/(1-y)) - Li2(-y(1-x)/(1-y))
    #     d/dx F = -log(1-x y)/(1-x) (the series derivative of Li_{1,1}); F(0,y)=0.
    #     Branch note: on 0<x,y<1 we have 0<1-x, 0<1-y, 0<(1-xy)/(1-y), and
    #     log((1-xy)/(1-y)) = log(1-xy) - log(1-y) (real log of positive args).
    z = -y*(1-x)/(1-y)
    dz = sp.diff(z, x)
    dF = (-sp.log(1-y)/(1-x)
          + (sp.log(1-x*y) - sp.log(1-y)) * dz / z)   # d/dx[-Li2(z)] = log(1-z) z'/z
    r = sp.factor(dF - (-sp.log(1-x*y)/(1-x)))
    print("  0.2  d/dx Li_11 closed form minus series derivative =", r)
    assert r == 0
    print("       F(0,y) = log(1)*log(1-y) + Li2(-y/(1-y)) - Li2(-y/(1-y)) = 0 (exact)")
    print("PART 0 OK")"""

assert old in content
content = content.replace(old, new)
with open(f'{wd}/verify_full.py','w') as f: f.write(content)
print("part0 patched")

# quick retest of part0 only
runner = 'import verify_full as V\nV.part0()\nprint("PART0 TRIAL OK")\n'
with open(f'{wd}/_trial0.py','w') as f: f.write(runner)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run([PY, f'{wd}/_trial0.py'], capture_output=True, text=True, cwd=wd, timeout=300)
print(r.stdout[-1500:])
print("STDERR:", r.stderr[-800:])
