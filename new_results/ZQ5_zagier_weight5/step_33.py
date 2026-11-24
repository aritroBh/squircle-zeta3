import sys, subprocess
s = open('steps/step_14_d8_lll.py').read()
# replace the whole comparison loop block with clean code
start = s.index('for tag, R in [("atoms"')
end = s.index('# use atom-level 830 relation as main')
newblock = '''for tag, R in [("atoms", (R1a, R2a)), ("kernel", (R1k, R2k))]:
    r1, r2 = R
    b1 = r1["rels"][0][1]; b2 = r2["rels"][0][1]
    m1 = r1["rels"][0][0]; m2 = r2["rels"][0][0]
    res1 = r1["rels"][0][2]; res2 = r2["rels"][0][2]
    p1 = prim(b1); p2 = prim(b2)
    same = p1 == p2
    print(tag, ": best rel max|a| P520=", m1, " P830=", m2, " identical after primitive: ", same)
    print("   resid P520=", mp.nstr(res1,2), " resid P830=", mp.nstr(res2,2))
    if not same:
        print("   P520:", p1)
        print("   P830:", p2)

'''
s = s[:start] + newblock + s[end:]
open('steps/step_14_d8_lll.py','w').write(s)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_14_d8_lll.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-2000:])
with open('progress.log','a') as f: f.write("step_14_d8_lll(r4): clean compare block; reran\n")
