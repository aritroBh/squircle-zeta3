import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
# fix checkpoint serialization: acb(str) may not parse str() output
src=open(os.path.join(m5b,'m5b_precision.py')).read()
src=src.replace("""    done = {}
    if part_path.exists():
        done = {k: acb(v) for k, v in json.loads(part_path.read_text()).items()}""",
"""    done = {}
    if part_path.exists():
        done = {k: acb(arb(v[0], v[1]), arb(v[2], v[3]))
                for k, v in json.loads(part_path.read_text()).items()}""")
src=src.replace("""        part_path.write_text(json.dumps({k: v.str(60) for k, v in done.items()}))""",
"""        part_path.write_text(json.dumps({k: [v.real.mid().str(70), v.real.rad().str(12),
                                             v.imag.mid().str(70), v.imag.rad().str(12)]
                                         for k, v in done.items()}))""")
open(os.path.join(m5b,'m5b_precision.py'),'w').write(src)
print('patched')
# quick check of arb/ acb roundtrip
from flint import acb, arb
x=acb(1.23456789)+acb(2).sqrt()*1e-3
s=[x.real.mid().str(70), x.real.rad().str(12), x.imag.mid().str(70), x.imag.rad().str(12)]
y=acb(arb(s[0],s[1]),arb(s[2],s[3]))
print(s)
print((x-y).str(10))
