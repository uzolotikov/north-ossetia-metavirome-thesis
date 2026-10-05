#!/usr/bin/env python3
"""Evaluate the two Biragzang.final match claims under explicit thresholds.

Inputs
  1. vOTU catalog vs Biragzang.final megablast (qseqid sseqid pident nident length
     qstart qend qlen slen bitscore)
  2. MV-4128 candidates vs Biragzang.final (12-column table of pipeline step 10)

For each query, non-overlapping HSPs are merged on the query axis so coverage is
not double counted, then the query is tested against a grid of thresholds.
"""
import sys
from collections import defaultdict


def merge(spans):
    out = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return sum(b - a + 1 for a, b in out)


def load(path, c_pid, c_aln, c_qs, c_qe, c_qlen, ncol):
    """-> {query: (qlen, [(pid, aln, qs, qe)])}"""
    d = defaultdict(list)
    qlen = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < ncol:
            continue
        q = p[0]
        qs, qe = sorted((int(p[c_qs]), int(p[c_qe])))
        d[q].append((float(p[c_pid]), int(p[c_aln]), qs, qe))
        qlen[q] = int(p[c_qlen])
    return {q: (qlen[q], v) for q, v in d.items()}


def report(label, data, n_total, grid):
    print(f"# {label}")
    print(f"#   queries with at least one HSP : {len(data)} of {n_total}")
    best_pid = max((h[0] for v in data.values() for h in v[1]), default=0)
    best_aln = max((h[1] for v in data.values() for h in v[1]), default=0)
    print(f"#   highest HSP identity anywhere : {best_pid:.2f} %")
    print(f"#   longest HSP anywhere          : {best_aln} bp")
    print("# min_pident\tmin_query_cov_pct\tqueries_passing")
    for pid, cov in grid:
        n = 0
        for q, (ql, hsps) in data.items():
            spans = [(a, b) for p_, al, a, b in hsps if p_ >= pid]
            if ql and 100.0 * merge(spans) / ql >= cov:
                n += 1
        print(f"# {pid}\t{cov}\t{n}")
    print("#")
    print("# per-query detail (best HSP by aligned length, any identity)")
    print("query\tqlen\tbest_pident\tbest_aln\tquery_cov_pct_any_identity")
    rows = []
    for q, (ql, hsps) in data.items():
        bp, ba, _, _ = max(hsps, key=lambda h: h[1])
        cov = 100.0 * merge([(a, b) for _, _, a, b in hsps]) / ql if ql else 0
        rows.append((ba, q, ql, bp, cov))
    for ba, q, ql, bp, cov in sorted(rows, reverse=True):
        print(f"{q}\t{ql}\t{bp:.2f}\t{ba}\t{cov:.2f}")
    print("#")


grid = [(95, 50), (95, 25), (95, 10), (95, 1), (90, 10), (80, 10), (75, 50)]

votu = load(sys.argv[1], 2, 4, 5, 6, 7, 10)
print("# Step 4b  do the vOTU catalog and the MV-4128 candidates match Biragzang.final?")
print("# Matches are counted for a grid of identity / query-coverage thresholds.")
print("#")
report("97 vOTU vs Biragzang.final (fresh megablast, unfiltered)", votu, 97, grid)

cand = load(sys.argv[2], 2, 3, 6, 7, 4, 12)
report("99 MV-4128 virus candidates vs Biragzang.final (pipeline step 10)", cand, 99, grid)
