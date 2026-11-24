import os, time
print("proof log:")
print(open('proof_div.log').read() if os.path.exists('proof_div.log') else "not yet created")
print("proof out:", open('run_proof_div.out').read() if os.path.exists('run_proof_div.out') else "none")
