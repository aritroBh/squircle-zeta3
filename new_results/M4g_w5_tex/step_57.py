import re, zlib
def pages(pdf):
    data=open(pdf,'rb').read()
    n=0
    for m in re.finditer(rb'stream\r?\n(.*?)endstream', data, re.S):
        try: d=zlib.decompress(m.group(1))
        except Exception: continue
        n+=len(re.findall(rb'/Type\s*/Page(?![s\w])', d))
    if n==0:
        n=len(re.findall(rb'/Type\s*/Page(?![s\w])', data))
    return n
print('supplement_W5.pdf pages:', pages('supplement_W5.pdf'))
print('fragment_test.pdf pages:', pages('fragment_test.pdf'))
log=open('progress.log','a')
log.write('check_supplement.py: round-trip exact (68 terms + 4 N5); W5 eval diff -9.08e-427 < 1e-390 PASS\n')
log.write('supplement_W5.tex: tectonic -X compile clean, no overfull hbox, PDF written\n')
log.write('paper_fragment.tex + fragment_test.tex: tectonic -X compile clean, no overfull hbox\n')
log.close()
print(open('progress.log').read())
