import subprocess
for pat in ["open('hexc2_common.py','w')", "open('numeric_cmzv.py','w')", "open('fast_cmzv.py','w')",
            "open('exact_linalg.py','w')", "open('verify_hexc2.py','w')", "open('result.md','w')",
            "member_modp", "fixed both"]:
    out = subprocess.run(['grep','-l','--',pat]+[f'step_{i:02d}.py' for i in range(1,75)],capture_output=True,text=True).stdout
    print(f'{pat}: {out.strip().splitlines()}')
print(open('progress.log').read())
