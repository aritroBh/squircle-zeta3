import os
os.chdir('<project>/research/08_theory/LEAD1_leading_genfun')
src = open('genfun.md').read()
old_start = src.index("### 3.5 Assembly")
old_end = src.index("## 4. Status")
new35 = """### 3.5 Assembly

For both parities the key-sum evaluates to `(d!/2) e_{d+1}`, so with d = w-1

    c(w) = delta_d/(d-1)! * (d!/2) e_{d+1} = delta_d * (d/2) * e_{d+1},

and since `delta_d = eps_{d-1} A_d/d!`, the prefactor is `delta_d d/2 = eps_{d-1} A_d/(2 (d-1)!)`.
With `d = w-1` this is exactly `c(w) = eps_{w-2} A_{w-1}/(2 (w-2)!) * e_w`, the stated law.

"""
src = src[:old_start] + new35 + src[old_end:]
open('genfun.md','w').write(src)
print('fixed 3.5')
