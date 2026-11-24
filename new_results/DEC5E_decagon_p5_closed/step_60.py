src = open('step_48_rfit.py').read()
# add D32 atoms to the setup and INV list
src = src.replace("""    atoms['K3ml2^2'] = K3m*l2**2; atoms['K3ml2lf'] = K3m*l2*lf; atoms['K3mlf^2'] = K3m*lf**2
    INV = [""", """    atoms['K3ml2^2'] = K3m*l2**2; atoms['K3ml2lf'] = K3m*l2*lf; atoms['K3mlf^2'] = K3m*lf**2
    li32 = pickle.load(open(os.path.join(here, 'li32_gold_800.pkl'),'rb'))
    atoms['D32a'] = mp.mpf(li32['vals'][0]) + lf**5/120
    atoms['D32b'] = mp.mpf(li32['vals'][1]) + lf**5/15
    atoms['D32c'] = mp.mpf(li32['vals'][2]) + mp.mpf(23)/120*lf**5
    atoms['D32d'] = mp.mpf(li32['vals'][3]) + mp.mpf(73)/120*lf**5
    INV = [""")
src = src.replace("'K3ml2^2','K3ml2lf','K3mlf^2']", "'K3ml2^2','K3ml2lf','K3mlf^2','D32a','D32b','D32c','D32d']")
src = src.replace("here = os.path.dirname(os.path.abspath(__file__))\n    W5", "here = os.path.dirname(os.path.abspath(__file__))\n    atoms0, INV0, lv0 = setup(800)  # ensure pickle path context\n    W5")
# fix: 'here' needed inside setup for li32 path -> define global HERE
src = src.replace("import mpmath as mp\nfrom q5_lll import fit, report",
                  "import mpmath as mp\nfrom q5_lll import fit, report\nHERE = os.path.dirname(os.path.abspath(__file__))")
src = src.replace("li32 = pickle.load(open(os.path.join(here, 'li32_gold_800.pkl'),'rb'))",
                  "li32 = pickle.load(open(os.path.join(HERE, 'li32_gold_800.pkl'),'rb'))")
src = src.replace("""    T1, T2 = int(sys.argv[3]), int(sys.argv[4])""", """    T1, T2 = int(sys.argv[3]), int(sys.argv[4])""")
open('step_48_rfit.py','w').write(src)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.time()
r = subprocess.run(['nice', PY, 'step_48_rfit.py', '450', '700', '520', '800'], capture_output=True, text=True, cwd='.')
print(r.stdout[-5000:])
if r.stderr: print("STDERR:", r.stderr[-900:])
print("wall %.0fs" % (time.time()-t0))
