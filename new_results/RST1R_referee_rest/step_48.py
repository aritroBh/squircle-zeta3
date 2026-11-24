code = open('step_09.py').read()
code = code.replace("""        Fe = mFH_exact(k, q, lam)
        Fg = FH_gauss(k, q, lamf); Fb = FH_box(k, q, lamf)
        rg = Fe/mpf(Fg.mid().str(40)); rb = Fe/mpf(Fb.mid().str(40))
        res[(k,tag)] = (Fe, rg, rb)
        print(f"  k={k:6d} {tag}: exact={mp.nstr(Fe,8)}  exact/gauss={mp.nstr(rg,8)}  exact/box={mp.nstr(rb,8)}", flush=True)""",
"""        Fe = mFH_exact(k, q, lam)
        Fg = FH_gauss(k, q, lamf); Fb = FH_box(k, q, lamf)
        Fe_f = arb(mp.nstr(Fe, 45))
        rg = (Fe_f/Fg); rb = (Fe_f/Fb)
        res[(k,tag)] = (Fe, rg, rb)
        print(f"  k={k:6d} {tag}: exact={mp.nstr(Fe,8)}  exact/gauss={rg.str(8)}  exact/box={rb.str(8)}", flush=True)""")
code = code.replace("""for kk,v in known.items():
    mine = res[kk][0]; ref = mpf(v)
    rel = abs(mine-ref)/ref
    print(f"  cross-check {kk}: mpmath={mp.nstr(mine,8)} arb={v} rel diff={mp.nstr(rel,3)}")
    assert rel < mpf('1e-6')""",
"""for kk,v in known.items():
    mine = res[kk][0]; ref = mpf(v)
    rel = abs(mine-ref)/ref
    print(f"  cross-check {kk}: mpmath={mp.nstr(mine,8)} arb={v} rel diff={mp.nstr(rel,3)}")
    assert rel < mpf('1e-6')""")
open('step_09.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_09.py'],capture_output=True,text=True,timeout=3600)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
