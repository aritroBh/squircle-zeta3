import subprocess, os
r = subprocess.run(['bash','-lc','ls ~/Library/TinyTeX/bin/ 2>/dev/null; ls ~/Library/TinyTeX/bin/*/pdflatex 2>/dev/null'], capture_output=True, text=True)
print(r.stdout, r.stderr)
WORK='<project>/research/08_theory/M4f_w5_extras'
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','gen_tex.py'],
                   cwd=WORK, capture_output=True, text=True, timeout=300)
print(r.stdout, r.stderr[-600:])
