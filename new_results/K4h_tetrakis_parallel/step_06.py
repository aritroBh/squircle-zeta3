"""step_06.py -- launch jobA_h768.py DETACHED (survives this session) at nice 10.

Command (exact, also in result.md):
  nohup nice -n 10 $VENV/python3 $TOOLS/detach.py $VENV/python3 \
      $CWD/jobA_h768.py > $CWD/jobA_h768_out.txt 2>&1 &
detach.py does os.setsid()+execvp, so the job gets its own session; the nice
level propagates across exec and is inherited by the 6 spawn workers.
"""
import subprocess, sys, os, time, json
CWD = "<project>/research/08_theory/K4h_tetrakis_parallel"
VENV = "<project>/.venv/bin/python3"
DETACH = "<project>/tools/detach.py"
JOB = os.path.join(CWD, "jobA_h768.py")

def main():
    assert os.path.exists(VENV), VENV
    assert os.path.exists(DETACH), DETACH
    assert os.path.exists(JOB), JOB
    out = open(os.path.join(CWD, "jobA_h768_out.txt"), "w")
    cmd = ["nohup", "nice", "-n", "10", VENV, DETACH, VENV, JOB]
    p = subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT,
                         stdin=subprocess.DEVNULL, cwd=CWD)
    print("launcher PID (detach.py, becomes the job after exec):", p.pid)
    print("command:", " ".join(cmd), "&")
    time.sleep(5)
    pidf = os.path.join(CWD, "jobA_h768.pid")
    if os.path.exists(pidf):
        print("job pid file:", open(pidf).read())
    json.dump({"launcher_pid": p.pid, "cmd": " ".join(cmd) + " &"},
              open(os.path.join(CWD, "launch_info.json"), "w"))

if __name__ == "__main__":
    main()
