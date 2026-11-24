#!/bin/bash
cd "<project>/research/08_theory/M5c_w6_high_precision"
nohup nice -n 10 "<project>/.venv/bin/python" -B m5b_precision.py w6 11500 11800 > w6_hp.out 2>&1 &
echo $! > w6_hp.pid
