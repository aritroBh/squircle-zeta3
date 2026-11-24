#!/bin/bash
cd "<project>/research/08_theory/M5c_w6_high_precision"
nohup nice -n 10 "<project>/.venv/bin/python" -B w6c_control.py 1300 > w6c_control.out 2>&1 &
echo $! > w6c_control.pid
