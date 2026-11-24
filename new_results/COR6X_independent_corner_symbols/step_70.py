import subprocess
for f in ['cor6_sym.py', 'cor6_pipe.py']:
    r = subprocess.run(["<project>/.venv/bin/python","-c",
        f"print(open('<project>/research/08_theory/COR6_true_point_symbols/{f}').read()[:5500])"],
        capture_output=True, text=True)
    print("=====", f); print(r.stdout[:5600])
