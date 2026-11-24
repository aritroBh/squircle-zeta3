import subprocess
for f in ["polylog_exact.py", "dilog_exact.py", "qsqrt2.py", "w2_integrals.py", "kummer.py",
          "<project>/research/08_theory/A3_W4_proof/w4_symbol.py"]:
    r = subprocess.run(["grep", "-n", "dps\|workdps\|prec", f], capture_output=True, text=True)
    print("=====", f)
    print(r.stdout)
