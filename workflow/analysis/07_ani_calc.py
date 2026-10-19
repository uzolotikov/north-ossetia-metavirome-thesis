#!/usr/bin/env python3
"""Нуклеотидная идентичность, взвешенная по длине, по таблице megablast: идентичность
MV-специфичного контаминанта между сборками MV-4125 и MV-4128.

Метод: HSP выбираются жадно по убыванию bitscore без перекрытий по оси запроса, так что каждое
основание запроса учитывается не более одного раза. Идентичность = sum(nident) / sum(длина выравнивания).

HSP фильтруются по идентичности MIN_PID и длине MIN_ALN. Фильтр важен: распределение идентичности HSP
двухвершинное — основная масса >=99,9 % и хвост ниже 97 % из повторов и паралогичных совпадений
с посторонними контигами MV-4125. Без фильтра хвост занижает оценку, и она перестаёт описывать
выравнивание ортологичных участков одного штамма. Основная оценка — при MIN_PID = 95, как принято
для ANI; таблица чувствительности к порогу печатается целиком, чтобы выбор можно было проверить.

Использование: 07_ani_calc.py <blast.tsv> <query.fna>
"""
import sys
from collections import defaultdict

MIN_PID = 95.0
MIN_ALN = 500

blast, query = sys.argv[1], sys.argv[2]

qlen, name = {}, None
with open(query) as fh:
    for line in fh:
        if line.startswith(">"):
            name = line[1:].split()[0]
            qlen[name] = 0
        else:
            qlen[name] += len(line.strip())

hsps = []
with open(blast) as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 14:
            continue
        qs, qe = sorted((int(p[7]), int(p[8])))
        hsps.append((p[0], float(p[13]), qs, qe, int(p[3]), int(p[4]), float(p[2])))


def greedy(min_pid, min_aln):
    byq = defaultdict(list)
    for q, bits, qs, qe, nid, aln, pid in hsps:
        if pid >= min_pid and aln >= min_aln:
            byq[q].append((bits, qs, qe, nid, aln))
    per, tot_aln, tot_id = {}, 0, 0
    for q, v in byq.items():
        taken, a_sum, i_sum = [], 0, 0
        for bits, qs, qe, nid, aln in sorted(v, reverse=True):
            if any(qs <= b and a <= qe for a, b in taken):
                continue
            taken.append((qs, qe))
            a_sum += aln
            i_sum += nid
        per[q] = (a_sum, i_sum)
        tot_aln += a_sum
        tot_id += i_sum
    return per, tot_aln, tot_id


print("# identity of the MV-specific contaminant between MV-4125 and MV-4128")
print("# query   : MV-4128 contigs >=50 kb (contam/ref_MV4128_ge50kb.fna)")
print("# subject : full MV-4125 metaSPAdes assembly")
print("# method  : unfiltered megablast; greedy non-overlapping HSPs on the query axis")
print("#")
print("# sensitivity to the HSP identity filter:")
print("# min_pid\tmin_aln\taligned_Mb\tidentity_pct\tdiff_per_100kb")
for mp in (0, 90, 95, 98, 99, 99.5):
    _, ta, ti = greedy(mp, MIN_ALN)
    idp = 100.0 * ti / ta
    print(f"# {mp}\t{MIN_ALN}\t{ta/1e6:.2f}\t{idp:.4f}\t{(100.0-idp)*1000:.1f}")
print("#")

per, tot_aln, tot_id = greedy(MIN_PID, MIN_ALN)
print(f"# HEADLINE at min_pid={MIN_PID}, min_aln={MIN_ALN}")
print("contig\tcontig_len\taligned_bp\tquery_cov_pct\tidentity_pct")
for q in sorted(qlen):
    a, i = per.get(q, (0, 0))
    pid = f"{100.0*i/a:.4f}" if a else "NA"
    print(f"{q}\t{qlen[q]}\t{a}\t{100.0*a/qlen[q]:.2f}\t{pid}")
ident = 100.0 * tot_id / tot_aln
print("#")
print(f"# contigs                : {len(qlen)}")
print(f"# query bases            : {sum(qlen.values())} ({sum(qlen.values())/1e6:.2f} Mb)")
print(f"# aligned bases          : {tot_aln} ({tot_aln/1e6:.2f} Mb)")
print(f"# identical bases        : {tot_id}")
print(f"# weighted identity      : {ident:.4f} %")
print(f"# differences per 100 kb : {(100.0-ident)*1000:.1f}")
