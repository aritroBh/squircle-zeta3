#!/usr/bin/env python3
"""Mechanical gate for one Lean target file.

Usage:  gate.py freeze <file.lean>     record the audited statement block
        gate.py check  <file.lean>     decide PROVED / FAIL, print JSON

File layout the gate enforces:

    import Mathlib
    ...
    -- BEGIN STATEMENTS
    namespace T01
    def stmt_1 : Prop := ...
    end T01
    -- END STATEMENTS
    theorem T01.proof_1 : T01.stmt_1 := by ...

`freeze` stores the text from the top of the file to the END marker (after an audit accepted it) in lean/audit/<name>.frozen.
`check` requires: frozen block unchanged; one `theorem <ns>.proof_k : <ns>.stmt_k` per `def stmt_k : Prop`;
no sorry/admit/native_decide/axiom/unsafe/implemented_by anywhere; the file compiles; and `#print axioms`
of every proof shows only propext, Classical.choice, Quot.sound.
"""
import json, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUDIT = os.path.join(ROOT, "audit")
OK_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}
FORBIDDEN = [r"\bsorry\b", r"\badmit\b", r"\bnative_decide\b", r"^\s*axiom\b", r"\bunsafe\b",
             r"implemented_by", r"\bopaque\b", r"set_option\s+debug", r"^\s*macro\b", r"^\s*elab\b",
             r"^\s*syntax\b", r"^\s*notation\b", r"run_cmd", r"#eval"]


def strip_comments(s):
    s = re.sub(r"/-.*?-/", "", s, flags=re.S)
    return re.sub(r"--.*", "", s)


def block(text):
    # everything from the top of the file to the END marker is frozen (imports and `open`s included)
    m = re.search(r"\A(.*?^-- BEGIN STATEMENTS\s*$.*?)^-- END STATEMENTS\s*$", text, flags=re.S | re.M)
    return m.group(1) if m else None


def frozen_path(path):
    return os.path.join(AUDIT, os.path.basename(path) + ".frozen")


def main():
    mode, path = sys.argv[1], os.path.abspath(sys.argv[2])
    text = open(path).read()
    b = block(text)
    out = {"file": os.path.relpath(path, ROOT), "ok": False, "problems": []}
    if b is None:
        out["problems"].append("no BEGIN/END STATEMENTS block")
        print(json.dumps(out, indent=1)); sys.exit(1)
    if mode == "freeze":
        os.makedirs(AUDIT, exist_ok=True)
        open(frozen_path(path), "w").write(b)
        out["ok"] = True
        print(json.dumps(out, indent=1)); return

    P = out["problems"]
    if not os.path.exists(frozen_path(path)):
        P.append("statement block was never frozen (no audit)")
    elif open(frozen_path(path)).read() != b:
        P.append("statement block differs from the audited one")
    code = strip_comments(text)
    for pat in FORBIDDEN:
        if re.search(pat, code, flags=re.M):
            P.append("forbidden: " + pat)
    ns = re.search(r"^namespace\s+(\S+)", strip_comments(b), flags=re.M)
    ns = ns.group(1) if ns else None
    if not ns:
        P.append("no namespace in statement block")
    stmts = re.findall(r"^\s*def\s+(stmt_\w+)\s*:\s*Prop\b", strip_comments(b), flags=re.M)
    if not stmts:
        P.append("no `def stmt_k : Prop` in statement block")
    after = strip_comments(text.split("-- END STATEMENTS", 1)[1])
    proofs = []
    for s in stmts:
        k = s[len("stmt_"):]
        if not re.search(r"^\s*theorem\s+%s\.proof_%s\s*:\s*%s\.%s\s*:=" % (re.escape(ns or ""), k,
                         re.escape(ns or ""), s), after, flags=re.M):
            P.append("missing `theorem %s.proof_%s : %s.%s :=`" % (ns, k, ns, s))
        else:
            proofs.append("%s.proof_%s" % (ns, k))
    out["statements"] = stmts
    if not P:
        with tempfile.NamedTemporaryFile("w", suffix=".lean", delete=False) as f:
            f.write(text + "\n" + "".join("#print axioms %s\n" % p for p in proofs))
            tmp = f.name
        try:
            r = subprocess.run(["lake", "env", "lean", tmp], cwd=ROOT, capture_output=True, text=True,
                               timeout=3600)
        finally:
            os.unlink(tmp)
        log = r.stdout + r.stderr
        if r.returncode != 0 or re.search(r"\berror\b", log):
            P.append("does not compile: " + log[-1500:])
        axioms = {}
        for p in proofs:
            m = re.search(r"'%s' depends on axioms: \[(.*?)\]" % re.escape(p), log, flags=re.S)
            if m:
                axioms[p] = [a.strip() for a in m.group(1).split(",")]
            elif re.search(r"'%s' does not depend on any axioms" % re.escape(p), log):
                axioms[p] = []
            else:
                P.append("no #print axioms output for " + p)
                continue
            bad = set(axioms[p]) - OK_AXIOMS
            if bad:
                P.append("%s uses axioms %s" % (p, sorted(bad)))
        out["axioms"] = axioms
    out["ok"] = not P
    print(json.dumps(out, indent=1))
    sys.exit(0 if out["ok"] else 1)


if __name__ == "__main__":
    main()
