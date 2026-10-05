#!/usr/bin/env python3
"""Производные таблицы каталога vOTU.

Использование: catalog_subsets.py <catalog_dir>
  representatives_ge10kb.fna       представители длиной >=10 000 б
  assembly_presence_ge5000.tsv     присутствие vOTU в сборках по членам кластера >=5 000 б
  assembly_presence_ge10000.tsv    то же по членам >=10 000 б (короткие члены присутствия не создают)
"""
import csv, pathlib, sys

P = pathlib.Path(sys.argv[1])
SAMPLES = ["MV-4122_S12", "MV-4125_S18", "MV-4128_S29"]

recs, name = [], None
for line in (P / "representatives_ge5kb.fna").read_text().splitlines(keepends=True):
    if line.startswith(">"):
        recs.append([line, []])
    else:
        recs[-1][1].append(line)
with (P / "representatives_ge10kb.fna").open("w") as f:
    for head, body in recs:
        if sum(len(l.strip()) for l in body) >= 10000:
            f.write(head + "".join(body))

rows = list(csv.DictReader((P / "membership.tsv").open(), delimiter="\t"))
order = list(dict.fromkeys(r["votu"] for r in rows))
for t in (5000, 10000):
    with (P / ("assembly_presence_ge%d.tsv" % t)).open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["votu"] + SAMPLES)
        for v in order:
            present = {r["sample"] for r in rows if r["votu"] == v and int(r["length"]) >= t}
            if present:
                w.writerow([v] + [int(s in present) for s in SAMPLES])
