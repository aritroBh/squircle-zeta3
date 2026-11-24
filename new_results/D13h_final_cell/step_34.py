import subprocess, sys
for f in ['step_18.py', 'step_19.py']:
    r = subprocess.run([sys.executable, f], capture_output=True, text=True, timeout=600)
    lines = r.stdout.splitlines()
    print(f"--- {f}: {len(lines)} lines ---")
    for l in lines[-4:]: print("  ", l)
    if r.stderr: print("  STDERR:", r.stderr[-200:])

# fix result.md references
r = open('result.md').read()
r = r.replace("29 cells, ratio <= 1.05 (`step_15.py`)", "29 cells, ratio <= 1.05 (`step_18.py`)")
r = r.replace("(`step_16.py`)", "(`step_19.py`)")
r = r.replace("  $PY step_15.py        # [100,400] ratio-1.05: 29 cells, sup 194.3809            (~1 s)",
              "  $PY step_18.py        # [100,400] ratio-1.05: 29 cells, sup 194.3809            (~1 s)")
r = r.replace("  $PY step_16.py        # spot-check planner sweep cells 0,20,41 (exact match)",
              "  $PY step_19.py        # spot-check planner sweep cells 0,20,41 (exact match)")
r = r.replace("step_15.py (tightening), step_16.py (sweep spot-check), step_17.py (junction + refinement),",
              "step_18.py (tightening), step_19.py (sweep spot-check), step_20.py (junction + refinement),")
open('result.md','w').write(r)
print("\nresult.md references fixed")
