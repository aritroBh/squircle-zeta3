import subprocess
print(subprocess.run(['bash','-c','head -c 400 "<project>/research/08_theory/E5d_nonlocal_pslq/b0_p5_convB_350.json"; echo; python3 -c "import json; d=json.load(open(\'<project>/research/08_theory/E5d_nonlocal_pslq/b0_p5_convB_350.json\')); print(type(d), list(d.keys())[:6] if isinstance(d,dict) else len(d))"'], capture_output=True, text=True).stdout)
