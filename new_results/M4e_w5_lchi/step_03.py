import os
for d in ['M4e_w5_lchi', 'M4c_w5_closed_form', 'M4d_w5_nice']:
    p = os.path.join('<project>/research/08_theory', d)
    print("="*70)
    print(d)
    for f in sorted(os.listdir(p)):
        fp = os.path.join(p, f)
        if os.path.isfile(fp):
            print("  ", f, os.path.getsize(fp))
        else:
            print("  DIR", f)
