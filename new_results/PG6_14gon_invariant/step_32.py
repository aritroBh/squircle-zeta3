base = '<project>/research/08_theory/'
pg = base+'PGEN1R_referee_alphabet/'
for i in [67,68,70,71,72,73,74,75,76]:
    fn = pg+f'step_{i}.py'
    try:
        c = open(fn).read()
        print("="*25, f"step_{i}", "="*25)
        print(c[:2400])
    except FileNotFoundError:
        pass
