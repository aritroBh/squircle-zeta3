#!/bin/bash
# resumable: verify.py skips targets that already have a RESULT line in verify.log
cd "$(dirname "$0")"
PY=../../../.venv/bin/python
ARGS=$($PY -c "from common import TARGETS; import sys; k=int(sys.argv[1]); print(' '.join(TARGETS[k::3]))" $1)
exec $PY verify.py $ARGS
