import sys, subprocess
lib = open('w5lib.py').read()
# add D=29 to letters_info
lib = lib.replace("""    else:
        raise ValueError(D)
    Sdesc={8:"S={2} (ramified; 2=(sqrt2)^2)",""",
"""    elif D==29:
        r=QD(0,1,29); eps=QD(Fraction(5,2),Fraction(1,2),29)
        letters=[eps,r,QD(2,0,29)]; names=["lE29","ls29","l2"]; disp=["log((5+sqrt29)/2)","log(sqrt29)","log(2)"]
        unit_idx=0; unit=eps; nonunit=[r,QD(2,0,29)]
        def pexps(N):
            if not _allowed(N,[2,29]): return None
            v2=_v_p(N,2)
            if v2%2!=0: return None
            return [_v_p(N,29), v2//2]
    else:
        raise ValueError(D)
    Sdesc={8:"S={2} (ramified; 2=(sqrt2)^2)",""")
lib = lib.replace("""           13:"S={13,2} (13=(sqrt13)^2 ramified, 2 inert)"}[D]""",
"""           13:"S={13,2} (13=(sqrt13)^2 ramified, 2 inert)",
           29:"S={29,2} (29=(sqrt29)^2 ramified, 2 inert)"}[D]""")
# chi for D=29 (generic Legendre for prime D congruent 1 mod 4)
lib = lib.replace("""    if D==13:
        if a%13==0: return 0
        r=pow(a%13,6,13); return 1 if r==1 else -1
    raise ValueError(D)""",
"""    if D in (13,29):
        if a%D==0: return 0
        r=pow(a%D,(D-1)//2,D); return 1 if r==1 else -1
    raise ValueError(D)""")
open('w5lib.py','w').write(lib)
p = open('pipe.py').read()
p = p.replace("    rad = {8:2, 5:5, 12:3, 13:13}[D]",
              "    rad = {8:2, 5:5, 12:3, 13:13, 29:29, 21:21, 17:17}[D]")
open('pipe.py','w').write(p)
src = r'''
import sys; sys.path.insert(0,'.')
import json, time
import mpmath as mp
import w5lib as W
import pipe
D=29
C, cands = pipe.phase_enum(D, B=8)
pipe.print_cands(D, cands, C)
K = pipe.phase_kernel(D, C, cands)
json.dump(dict(K=K), open(pipe.cand_dir(D)+"/kernel.json","w"))
t0=time.time()
vals = pipe.phase_numerics(D, C, cands, P=1150)
print(f"numerics {time.time()-t0:.0f}s")
for P in (520, 830):
    R = pipe.phase_lll(D, C, cands, K, vals, P=P, use_atoms=False)
    print(f"P={P}: exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])}")
    for mx,a,res in R["exact_rels"][:2]:
        print("   a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
'''
open('steps/step_46_d29.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_46_d29.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print("STDERR:", r.stderr[-1500:])
with open('progress.log','a') as f: f.write("step_46_d29: D=29 pipeline through kernel LLL\n")
