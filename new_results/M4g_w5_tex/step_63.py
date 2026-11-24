import subprocess, re, zlib
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'check_supplement.py','../M4f_w5_extras/R5_best.json','supplement_W5_best.tex'],capture_output=True,text=True,timeout=1500)
print(r.stdout[-500:]); print(r.stderr[-400:])
r=subprocess.run(['nice','-n','10','/opt/homebrew/bin/tectonic','-X','compile','supplement_W5_best.tex'],capture_output=True,text=True,timeout=600)
print('--- tectonic best ---'); print(r.stdout[-400:]); print(r.stderr[-400:])
def pages(pdf):
    data=open(pdf,'rb').read(); n=0
    for m in re.finditer(rb'stream\r?\n(.*?)endstream', data, re.S):
        try: d=zlib.decompress(m.group(1))
        except Exception: continue
        n+=len(re.findall(rb'/Type\s*/Page(?![s\w])', d))
    return n or len(re.findall(rb'/Type\s*/Page(?![s\w])', data))
for f in ('supplement_W5.pdf','supplement_W5_best.pdf','fragment_test.pdf'):
    print(f,'pages:',pages(f))
log=open('progress.log','a')
log.write('gen_supplement.py generic (argv, dynamic groups, +4 ARG_FORMS, first-line + stripped); diff vs old = 3 sign-only lines\n')
log.write('check_supplement.py generic: R5_nice_reduced round-trip exact, W5 diff -9.08e-427 PASS\n')
log.write('R5_best.json: gen+check ran on supplement_W5_best.tex; tectonic compiled both supplements clean\n')
log.close()
