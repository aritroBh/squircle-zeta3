code = open('step_12.py').read()
code = code.replace("""from flint import arb
PI = arb.pi()""",
"""from flint import arb
def mpf_of(x): return mp.mpf(str(x.mid().str(30)).split(' ')[0].lstrip('['))
PI = arb.pi()""")
code = code.replace("mp.mpf(S_exact)/mp.mpf(S_bound.mid().str(30))", "mp.mpf(S_exact)/mpf_of(S_bound)")
code = code.replace("assert S_exact < mp.mpf(S_bound.mid().str(30))", "assert S_exact < mpf_of(S_bound)")
code = code.replace("""print(f"    q2 : exact={mp.nstr(exact_terms[0]+2*exact_terms[1],10)} bound={(b_q2/(8/PI)).str(10)}")""",
"""print(f"    q2 : exact={mp.nstr(exact_terms[0]+2*exact_terms[1],10)} bound={(b_q2/(8/PI)).str(10)}")""")
code = code.replace("mp.nstr(mp.mpf(T2e)/mp.mpf(T2b.mid().str(30)),6)", "mp.nstr(mp.mpf(T2e)/mpf_of(T2b),6)")
code = code.replace("assert T2e < mp.mpf(T2b.mid().str(30))", "assert T2e < mpf_of(T2b)")
open('step_12.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'step_12.py'],capture_output=True,text=True,timeout=3600)
print("rc:",r.returncode); print(r.stdout); print(r.stderr[-1500:])
