for n in (23, 32, 31):
    p = f'<project>/research/08_theory/M4b_w4pp_exact/step_{n}.py'
    try:
        c = open(p).read()
        print('='*20, n, len(c))
        print(repr(c[:200]))
    except FileNotFoundError:
        print(n, 'missing')
