import os
steps = sorted((f for f in os.listdir('.') if f.startswith('step_')), key=lambda s: int(s.split('_')[1].split('.')[0]))
for f in steps:
    src = open(f).read()
    first = [l for l in src.split('\n') if l.strip() and not l.startswith('import')][:2]
    # heuristics
    tag = ''
    for key, t in [('enum_maps13()', 'N=6 maps'), ('def map_count', 'map counts'), ('cross_values13', 'cross-ratios'),
                   ('make_relations(n, lower_t)', 'relsets'), ('residual_fmpq', 'cert exact'), ('tag_families', 'rows numeric/families'),
                   ('my13_hpl_to_cmzv', 'conversion'), ('seg_wrong', 'wrong order'), ('hpl_wrong', 'wrong order'),
                   ('letterform_to_poly', '1-forms exact'), ('id2 = R32', 'ID 250dps'), ('M3_closed', 'M3 1150'),
                   ('rank_modp_fast', 'dims'), ('verify_cor13c', 'their verifier'), ('y_bad', 'tamper'),
                   ('verify_cor13cr', 'my verifier'), ('disp_w3', 'displays'), ('reflib13.py', 'lib written'),
                   ('mymaps.py', 'mymaps written'), ('cross_reps13(Fr(z))', 'frames'), ('T4_hcf', 'frames'),
                   ('my_target2(mi)', 'map independence')]:
        if key in src: tag = t
    print(f, '|', tag, '|', (first[0][:70] if first else ''))
