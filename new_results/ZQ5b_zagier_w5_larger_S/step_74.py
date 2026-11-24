import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
with open('progress.log','a') as f:
    f.write('step_08_d13_enum: 487 exceptional S-units for D=13 S={13,2,3} (B=8) saved out/d13/\n')
    f.write('step_09g_d13_kernel: weight-5 kernel dim 263 (rank 224, N=487), unsaturated (saturation too slow, not needed for search); LLL-reduced; companions saved out/d13/kernel5.json\n')
print("logged")
