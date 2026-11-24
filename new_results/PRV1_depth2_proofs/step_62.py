wd = 'research/08_theory/PRV1_depth2_proofs'
with open(f'{wd}/verify_full.py') as f:
    content = f.read()

chain = r'''

# ======================================================================
# Chain-rule assembly (exact): the dG, dH, dGd, dHd formulas used in
# Parts 2 and 3 are derived HERE from Part 0's partial-derivative formulas
# and the Li_{1,1} closed form, by exact rational simplification on the curves.
# ======================================================================
def part_chainrule():
    print("=" * 74)
    print("CHAIN-RULE ASSEMBLY (exact): dG, dH, dGd, dHd from Part 0")
    tt = sp.symbols('t', positive=True)
    def li11(xe, ye, log1x, log1y, L2a, L2b):
        # Li11(xe,ye) via Part 0.2, with the exact argument evaluations supplied
        return log1x*log1y + L2a - L2b
    def dLi21_y(xe, ye, L2x, L2xy):
        return (L2x - L2xy/ye)/(1-ye)
    # ---- octagon G: (x,y) = (-t, -(1-t)/(2t)) ----
    xe, ye = -tt, -(1-tt)/(2*tt)
    # exact argument identities (rational, checked):
    ids = {"1-x = 1+t": sp.simplify(1-xe - (1+tt)),
           "1-y = (1+t)/(2t)": sp.simplify(1-ye - (1+tt)/(2*tt)),
           "-y/(1-y) = (1-t)/(1+t) = c": sp.simplify(-ye/(1-ye) - (1-tt)/(1+tt)),
           "-y(1-x)/(1-y) = 1-t": sp.simplify(-ye*(1-xe)/(1-ye) - (1-tt)),
           "xy = (1-t)/2 = b": sp.simplify(xe*ye - (1-tt)/2)}
    # ---- octagon H: (x,y) = ((1-t)/2, -2t/(1-t)) ----
    xeH, yeH = (1-tt)/2, -2*tt/(1-tt)
    ids.update({"H 1-x = (1+t)/2 = h": sp.simplify(1-xeH - (1+tt)/2),
                "H 1-y = (1+t)/(1-t)": sp.simplify(1-yeH - (1+tt)/(1-tt)),
                "H -y/(1-y) = 2t/(1+t) = u": sp.simplify(-yeH/(1-yeH) - 2*tt/(1+tt)),
                "H -y(1-x)/(1-y) = t": sp.simplify(-yeH*(1-xeH)/(1-yeH) - tt),
                "H xy = -t": sp.simplify(xeH*yeH + tt)})
    # ---- decagon G: (x,y) = (-t, -(1-t)/t) ----
    xeG2, yeG2 = -tt, -(1-tt)/tt
    ids.update({"dG 1-x = 1+t": sp.simplify(1-xeG2 - (1+tt)),
                "dG 1-y = 1/t": sp.simplify(1-yeG2 - 1/tt),
                "dG -y/(1-y) = 1-t": sp.simplify(-yeG2/(1-yeG2) - (1-tt)),
                "dG -y(1-x)/(1-y) = 1-t^2": sp.simplify(-yeG2*(1-xeG2)/(1-yeG2) - (1-tt**2)),
                "dG xy = 1-t": sp.simplify(xeG2*yeG2 - (1-tt))})
    # ---- decagon H: (x,y) = (1-t, -t/(1-t)) ----
    xeH2, yeH2 = 1-tt, -tt/(1-tt)
    ids.update({"dH 1-x = t": sp.simplify(1-xeH2 - tt),
                "dH 1-y = 1/(1-t)": sp.simplify(1-yeH2 - 1/(1-tt)),
                "dH -y/(1-y) = t": sp.simplify(-yeH2/(1-yeH2) - tt),
                "dH -y(1-x)/(1-y) = t^2": sp.simplify(-yeH2*(1-xeH2)/(1-yeH2) - tt**2),
                "dH xy = -t": sp.simplify(xeH2*yeH2 + tt)})
    for name, r in ids.items():
        print("   ", name, ": residual =", r)
        assert r == 0
    # assemble derivatives (real-log splits of positive arguments on 0<t<1):
    # octagon G: log(1-x)=B, log(1-y)=B-C-R=q-R, L2(-y/(1-y))=Dc, L2(-y(1-x)/(1-y))=D1
    Li11_G = B*(q-R) + Dc - D1
    dG_asm = (Li11_G/xe)*sp.diff(xe, tt) + dLi21_y(xe, ye, Dm, Db)*sp.diff(ye, tt)
    r = sp.factor(dG_asm - dG)
    print("    assembled dG - stated dG =", r); assert r == 0
    # octagon H: log(1-x)=q, log(1-y)=-lc, L2(u)=Du, L2(t)=Dt
    Li11_H = -lc*q + Du - Dt
    dH_asm = (Li11_H/xeH)*sp.diff(xeH, tt) + dLi21_y(xeH, yeH, Db, Dm)*sp.diff(yeH, tt)
    r = sp.factor(dH_asm - dH)
    print("    assembled dH - stated dH =", r); assert r == 0
    # decagon G: log(1-x)=B, log(1-y)=-R, L2(1-t)=D1, L2(1-t^2)=D1t2
    Li11_Gd = B*(-R) + D1 - D1t2
    dGd_asm = (Li11_Gd/xeG2)*sp.diff(xeG2, tt) + dLi21_y(xeG2, yeG2, Dm, D1)*sp.diff(yeG2, tt)
    r = sp.factor(dGd_asm - dGd)
    print("    assembled dGd - stated dGd =", r); assert r == 0
    # decagon H: log(1-x)=R, log(1-y)=-A, L2(t)=Dt, L2(t^2)=Dt2
    Li11_Hd = R*(-A) + Dt - Dt2
    dHd_asm = (Li11_Hd/xeH2)*sp.diff(xeH2, tt) + dLi21_y(xeH2, yeH2, D1, Dm)*sp.diff(yeH2, tt)
    r = sp.factor(dHd_asm - dHd)
    print("    assembled dHd - stated dHd =", r); assert r == 0
    print("CHAIN-RULE ASSEMBLY OK")
'''

# insert the chain function before "def main()" and add the call
old_main = """def main():
    part0()
    part1()
    part2()"""
new_main = """def main():
    part0()
    part1()
    part_chainrule()
    part2()"""
assert old_main in content
content = content.replace(old_main, new_main)
# insert chain function definition just before part4 section header
anchor = "\n\n# ======================================================================\n# Part 4: independent numeric evidence"
assert anchor in content
content = content.replace(anchor, chain + anchor)
with open(f'{wd}/verify_full.py','w') as f: f.write(content)
print("chain-rule assembly inserted")

import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run([PY, f'{wd}/verify_full.py'], capture_output=True, text=True, cwd=wd, timeout=1800)
out = r.stdout
i = out.find('CHAIN-RULE')
print(out[i-80:i+2600])
print("TAIL:", out[-120:])
print("STDERR:", r.stderr[-1200:])
with open(f'{wd}/verify_full.out','w') as f:
    f.write(out + "\nSTDERR:\n" + r.stderr)
