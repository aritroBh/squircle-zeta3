import os
os.chdir("<project>")
for root, dirs, files in os.walk("."):
    dirs[:] = [d for d in dirs if d not in ('.venv', '.git', '__pycache__')]
    for f in files:
        if f in ('prev_result.md', 'esum_engine.py', 'Jtab_t.pkl', 'E3_general.pkl', 'atom_values.json'):
            print(os.path.join(root, f))
