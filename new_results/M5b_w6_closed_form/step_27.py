import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
src=open(os.path.join(m5b,'m5b_precision.py')).read()
old = """    part_path = HERE/f'{target}_{bits}_partial.json'
    done = {}
    if part_path.exists():
        done = {k: acb(arb(v[0], v[1]), arb(v[2], v[3]))
                for k, v in json.loads(part_path.read_text()).items()}"""
new = """    part_path = HERE/f'{target}_{bits}_partial.json'
    full_digits = int(bits*0.301029995663981)+10
    def pack(v):
        return [v.real.mid().str(full_digits, radius=False), str(v.real.rad().upper()),
                v.imag.mid().str(full_digits, radius=False), str(v.imag.rad().upper())]
    def unpack(v):
        return acb(arb(v[0], arb(v[1])), arb(v[2], arb(v[3])))
    done = {}
    if part_path.exists():
        done = {k: unpack(v) for k, v in json.loads(part_path.read_text()).items()}"""
assert old in src
src = src.replace(old, new)
old2 = """        part_path.write_text(json.dumps({k: [v.real.mid().str(70), v.real.rad().str(12),
                                             v.imag.mid().str(70), v.imag.rad().str(12)]
                                         for k, v in done.items()}))"""
new2 = """        part_path.write_text(json.dumps({k: pack(v) for k, v in done.items()}))"""
assert old2 in src
src = src.replace(old2, new2)
open(os.path.join(m5b,'m5b_precision.py'),'w').write(src)
# test roundtrip
from flint import acb, arb
x = acb(arb('1.23456789','1e-25'), arb('3e-30'))
packed = [x.real.mid().str(200, radius=False), str(x.real.rad().upper()),
          x.imag.mid().str(200, radius=False), str(x.imag.rad().upper())]
y = acb(arb(packed[0], arb(packed[1])), arb(packed[2], arb(packed[3])))
print(packed)
print('contains orig diff:', (x-y).abs() < 3e-25, (x-y).str(5))
print('OK')
