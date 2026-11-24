import os
open('progress.log','a').write('DONE: report composed; all tasks 1-6 complete\n')
print(open('progress.log').read())
print(sorted(os.listdir('steps')))
