import os
os.chdir('<project>/research/08_theory/DEC5E_decagon_p5_closed')
build = '''"""step_28_build_tables.py - build decagon G-tables at (dps,N) given by argv, with self-checks."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gengine_dec2p
import mpmath as mp

def main():
    dps, N = int(sys.argv[1]), int(sys.argv[2])
    procs = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    t0 = time.time()
    G1, letters = gengine_dec2p.gtab(dps, N, cachedir=os.path.dirname(os.path.abspath(__file__)), procs=procs)
    mp.mp.dps = dps
    el = time.time()-t0
    L0, LP, LQ = letters['L0'], letters['LP'], letters['LQ']
    r1 = abs(G1[(L0,)*4+(LP,)] + mp.polylog(5, 1/LP))
    r2 = abs(G1[(L0,)*4+(LQ,)] + mp.polylog(5, 1/LQ))
    sh = abs(G1[(LP,)]*G1[(LQ,)] - G1[(LP,LQ)] - G1[(LQ,LP)])
    print("built dps=%d N=%d procs=%d in %.1fs, %d words" % (dps, N, procs, el, len(G1)))
    print("G(0^4,-phi;1)+Li5(-1/phi):", mp.nstr(r1, 3))
    print("G(0^4,phi^2;1)+Li5(1/phi^2):", mp.nstr(r2, 3))
    print("shuffle G(p)G(q):", mp.nstr(sh, 3))

if __name__ == "__main__":
    main()
'''
open('step_28_build_tables.py','w').write(build)
# fix gengine_dec2p.py so spawn workers import it fine (it's a module import, safe)
import subprocess
PY = '<project>/.venv/bin/python'
t0=__import__('time').time()
r = subprocess.run(['nice', PY, 'step_28_build_tables.py', '120', '380', '10'], capture_output=True, text=True, cwd='.')
print(r.stdout, r.stderr[-500:] if r.stderr else "")
print("wall %.1fs" % (__import__('time').time()-t0))
