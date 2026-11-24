#!/bin/zsh
cd "<project>/research/08_theory/CX7_tow1_review"
echo "launch $(date)" >> codex.log
codex exec --skip-git-repo-check -C "$PWD" -s workspace-write -m gpt-6.1-sol -c model_reasoning_effort=high -o "$PWD/codex_final.md" - < brief.md >> codex.log 2>&1
echo "exit=$? $(date)" >> codex.log
