for f in ['step_01.py','step_02.py','step_03.py','step_04.py','step_05.py']:
    print("="*30, f, "="*30)
    with open(f) as fh:
        print(fh.read())
