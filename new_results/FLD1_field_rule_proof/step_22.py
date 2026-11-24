WD='<project>/research/08_theory/FLD1_field_rule_proof'
s=open(WD+'/step_21_fastmu.py').read()
# Replace the verification block with an arity-normalized one
old_start = s.index("    # verify k<=3 against G1 table")
old_end = s.index("    json.dump(")
newblock = '''    # verify k<=3 against G1 table flat_exact_N320_K3.json (n<=100); normalize arity
    G1='<project>/research/10_frontier/G1_zeta3_any_angle/'
    ref = json.load(open(G1+'flat_exact_N320_K3.json'))
    def normkey(kk):
        comps = [int(x) for x in kk.strip('()').split(',')]
        comps = (comps + [0,0,0,0])[:4]
        return ','.join(map(str,comps))
    ok=True
    for k in (0,1,2,3):
        for n in range(1,101):
            mine = table[k][n]
            theirs = {normkey(kk): vv for kk,vv in ref[str(k)][str(n)].items()}
            for kk in set(list(mine)+list(theirs)):
                if Fr(mine.get(kk,'0')) != Fr(theirs.get(kk,'0')):
                    ok=False; print('MISMATCH',k,n,kk,mine.get(kk),theirs.get(kk))
    print('k<=3 matches G1 flat_exact table (n<=100):', ok)
'''
s = s[:old_start] + newblock + s[old_end:]
open(WD+'/step_21_fastmu.py','w').write(s)
import subprocess
VENV='<project>/.venv/bin/python'
r = subprocess.run(['nice',VENV,'step_21_fastmu.py','5','100'],capture_output=True,text=True,cwd=WD)
for line in r.stdout.splitlines():
    if 'MISMATCH' in line: continue
    print(line)
print(r.stderr[-1000:] if r.returncode else '')