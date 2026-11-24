import datetime
wd = '<project>/research/08_theory/TOW5R_review'
with open(wd+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now():%H:%M:%S} step_24: result.md final written; TOW5R review complete\n")
print('logged')
