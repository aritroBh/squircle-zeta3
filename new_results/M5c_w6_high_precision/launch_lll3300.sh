#!/bin/bash
cd "<project>/research/08_theory/M5c_w6_high_precision"
nohup nice -n 10 "<project>/.venv/bin/python" -B w6c_lll.py 3300 > w6c_lll_3300.out 2>&1 &
echo $! > w6c_lll_3300.pid
