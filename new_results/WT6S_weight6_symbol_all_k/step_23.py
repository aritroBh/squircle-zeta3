
import os, shutil, datetime
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
WD = os.path.join(T, 'WT6S_weight6_symbol_all_k')
os.chdir(WD)
print("CWD:", os.getcwd())

# Copy modules we reuse (never edit elsewhere)
copies = [
    (os.path.join(ROOT, 'paper/squircle_note/scripts/qsqrt2.py'), 'qsqrt2.py'),
    (os.path.join(T, 'A3_W4_proof/w4_symbol.py'), 'w4_symbol.py'),
    (os.path.join(T, 'TOWN1_tower_large_n/claude_checks/deg5_residue_polynomials.json'), 'deg5_residue_polynomials.json'),
    (os.path.join(T, 'TOWN1_tower_large_n/claude_checks/deg4_residue_polynomials.json'), 'deg4_residue_polynomials.json'),
    (os.path.join(T, 'WTS_weight_k_structure/theory/kern_5.pkl'), 'kern_5.pkl'),
    (os.path.join(T, 'WT5_weight5_all_k/kernels_deg4.json'), 'kernels_deg4.json'),
    (os.path.join(T, 'WT5_weight5_all_k/deltas_deg4.json'), 'deltas_deg4.json'),
    (os.path.join(T, 'M5a_next_coefficient/omega5.json'), 'omega5.json'),
]
for src, dst in copies:
    shutil.copy(src, os.path.join(WD, dst))
    print("copied", dst)

# Start progress log
with open('progress.log', 'a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} setup: copied qsqrt2.py, w4_symbol.py, deg5/deg4 residue polynomials, kern_5.pkl, kernels_deg4.json, deltas_deg4.json, omega5.json\n")
print("done")
