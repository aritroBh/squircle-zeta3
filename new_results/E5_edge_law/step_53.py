code = open('steps/step_09.py').read()
code = code.replace("# kinks: j0 = 0 (n=(0,1) through V0), j1 = 1 (through V1); non-local j = 2,3,4",
                    "# kinks: j0 = 0 (n=(0,1) through V0), j1 = m-1 = 4 (through V1); non-local j = 1,2,3")
code = code.replace("alpha_ = aj(0, 1); beta_ = aj(1, 0)", "alpha_ = aj(0, 1); beta_ = aj(m-1, 0)")
code = code.replace("           l0s: lj[0], l1s: lj[1]}", "           l0s: lj[0], l1s: lj[m-1]}")
code = code.replace("            sub[taus[k-1]] = mp.fsum(aj(j, t)*lj[j]**k for j in range(2, m))\n            sub[gams[k-1]] = mp.fsum(ap[j]*lj[j]**k for j in range(2, m))",
                    "            sub[taus[k-1]] = mp.fsum(aj(j, t)*lj[j]**k for j in range(1, m-1))\n            sub[gams[k-1]] = mp.fsum(ap[j]*lj[j]**k for j in range(1, m-1))")
open('steps/step_09.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_09.py'],capture_output=True,text=True,timeout=3000)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-2500:])
