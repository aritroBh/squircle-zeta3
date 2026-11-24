base='<project>'
for f in ['research/08_theory/A3h_weight4_target/result.md',
          'research/08_theory/A3h_weight4_target/progress.log']:
    print('#####',f)
    print(open(base+'/'+f).read())
    print()
