import subprocess, os
PY='<project>/.venv/bin/python'
base='<project>'
# List key directories
for d in ['research/08_theory/A3h_weight4_target',
          'research/08_theory/A3c_claude_elementary/Q',
          'research/08_theory/A3f_level8_relations',
          'research/08_theory/A3f4_glanois_basis']:
    p=os.path.join(base,d)
    print('===',d)
    try:
        print(subprocess.run(['ls','-la',p],capture_output=True,text=True).stdout)
    except Exception as e:
        print(e)
