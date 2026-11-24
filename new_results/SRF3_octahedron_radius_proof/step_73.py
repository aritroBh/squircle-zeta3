
import os, datetime
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
print("result.md size:", os.path.getsize(wd+"/result.md"))  # check it was written
log = open(wd+"/progress.log","a")
msg = "[{}] SRF3R referee repairs written into proof.md + result.md: V-split (3.2), direct non-cancellation (3.3, H0=sqrt(1+2e^(2w))), d-display sign -sqrt(3pi/4), tau* sign -, phase typos, open-items repair paths; d-sign verified: ratios about -1.53499\n".format(datetime.date.today())
log.write(msg)
log.close()
print(open(wd+"/progress.log").read()[-300:])
