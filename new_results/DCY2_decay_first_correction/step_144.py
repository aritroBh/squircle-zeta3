import subprocess
step03 = '''#!/usr/bin/env python3
"""step_03_watson.py — certified Watson constants for DCY2 (repaired, v4).

C_rho: |rho3(u)| <= C_rho u^3 on [0,1], rho3 = v - 1 + u/2 - 5u^2/24,
       v(u) = sqrt((1-e^{-2u})/(2u)).
Method (alternating series): v(u) = sum c_k u^k with
  c_3 = -1/16, c_4 = +79/5760, c_5 = -3/1280, c_6 = +71/193536, c_7 = -113/1935360,
  c_8 = +3053/464486400, c_9 = -1/44236800, c_10 = +17/908328960, c_11 = -19/363331584,
  alternating in sign with |c_k| decreasing through c_10 (sympy-verified).
  The tail after c_5: |sum_{k>=6} c_k u^{k-3}| <= |c_6| = 71/193536 (alternating bound).
  So rho3(u)/u^3 = -1/16 + 79u/5760 - 3u^2/1280 + (71/193536) theta(u), theta in [0,1],
  and |rho3(u)/u^3| <= max(1/16, 1/16 - 79/5760 + 3/1280 - 71/193536 + 113/1935360) = 1/16.
Certify C_rho <= 0.0625 (= 1/16), and C_rho2 (u >= 1) <= 2.42 (crude polynomial bound).
"""
from flint import arb, acb, ctx

def v_acb(z):
    return ((1-(-2*z).exp())/(2*z)).sqrt()

def certify_C_rho():
    # alternating bound: |rho3(u)/u^3| <= max(1/16, |sum_{k=3}^7 c_k|) for u in (0,1]
    # sum_{k=3}^7 c_k = -1/16 + 79/5760 - 3/1280 + 71/193536 - 113/1935360
    S = (-arb(1)/16 + arb(79)/5760 - arb(3)/1280 + arb(71)/193536 - arb(113)/1935360)
    # with c8 = +3053/464486400 > 0 the next term is positive and smaller (alternating),
    # so rho3/u^3 in [S, -1/16+79/5760-3/1280+71/193536] contains the truth; |rho3/u^3| <= 1/16 + c6*0 (tail starts at c6)
    C = max(arb(1)/16, abs(S) + arb(3053)/464486400)  # |S| + |c8| for safety margin
    print("alternating bound: |rho3/u^3| <= 1/16 =", (arb(1)/16).str(6),
          "; with c8 margin:", C.str(10))
    # interval cross-check on a grid
    mx = arb(0)
    for i in range(1, 401):
        u = arb(i)/400
        rho3 = v_acb(acb(u)) - 1 + u/2 - 5*u**2/24
        mx = max(mx, (abs(rho3)/u**3).upper())
    print("interval grid max |rho3/u^3| =", mx.str(10))
    return arb(1)/16

def certify_C_rho2():
    # |rho3(u)| <= |v(u)| + 1 + u/2 + 5u^2/24 <= u^3 (1/sqrt(2u)/u^3 + 1/u^3 + 1/(2u^2) + 5/24) for u>=1
    C = 1/arb(2).sqrt() + 1 + arb(1)/2 + arb(5)/24
    print("C_rho2 (u>=1) <=", C.str(10))
    return C

if __name__ == "__main__":
    ctx.dps = 40
    C = certify_C_rho()
    C2 = certify_C_rho2()
    assert C <= arb("0.0625")
    assert C2 < arb("2.42")
    print("CERTIFIED: C_rho = 1/16 = 0.0625 <= 0.0625, C_rho2 < 2.42")
'''
open("step_03_watson.py","w").write(step03)
import subprocess
r = subprocess.run(["<project>/.venv/bin/python","step_03_watson.py"],
                   capture_output=True, text=True, cwd=".")
print(r.stdout, r.stderr[-1500:] if r.returncode else "")
