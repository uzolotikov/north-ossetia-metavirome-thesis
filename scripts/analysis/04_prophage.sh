#!/bin/bash
set -uo pipefail
K=/NatureUsers/nsinyavskiy/metaviroms_kgt; W=/NatureUsers/glitch/mv_review
B=/NatureUsers/nsinyavskiy/miniconda3/envs/blast-env/bin; M=/NatureUsers/nsinyavskiy/miniconda3/envs/mapping-env/bin
cd $W/prophage
$B/makeblastdb -in $W/step6/ref_derep.fna -dbtype nucl -out ref_derep > /dev/null
$B/blastn -task megablast -query $K/votu/no_NC/catalog/votu_no_NC_95_85/representatives_ge5kb.fna -db ref_derep -perc_identity 99 -evalue 1e-50 -num_threads 16 -outfmt "6 qseqid sseqid pident length qstart qend sstart send qlen slen" > votu_vs_ref.tsv
python3 - <<'P'
import collections
best={}
hs=collections.defaultdict(list)
for l in open("votu_vs_ref.tsv"):
    q,s,p,L,qs,qe,ss,se,ql,sl=l.split("\t"); hs[(q,s)].append((int(ss),int(se),int(L),int(ql),int(sl)))
for (q,s),h in hs.items():
    al=sum(x[2] for x in h); ql=h[0][3]; sl=h[0][4]
    if al/ql<0.9 or sl<ql+5000: continue
    lo=min(min(x[0],x[1]) for x in h); hi=max(max(x[0],x[1]) for x in h)
    if q not in best or al>best[q][3]: best[q]=(s,lo,hi,al,sl)
with open("embedded.tsv","w") as o:
    for q,(s,lo,hi,al,sl) in sorted(best.items()): o.write(f"{q}\t{s}\t{lo}\t{hi}\t{sl}\n")
print("embedded vOTUs:",len(best))
P
LIBS="MV-4122_S12 MV-4125_S18 MV-4128_S29 4122-0-22_S11 4122-0-44_S10 4125-0-22_S17 4125-0-44_S16 4128-0-22_S28 4128-0-44_S27 4128_S30 NC1_S64"
echo -e "votu\thost_contig\tstart\tend\thost_len\tlibrary\tmedian_depth_in\tmedian_depth_flanks\tratio_in_over_flank" > prophage_ratios.tsv
while IFS=$'\t' read v c lo hi sl; do
  for lib in $LIBS; do
    $M/samtools depth -a -r "$c" $W/step6/bam/$lib.bam 2>/dev/null | python3 -c "
import sys,statistics
lo,hi=$lo,$hi; pad=1000; ins=[];fl=[]
for l in sys.stdin:
    f=l.split('\t'); p=int(f[1]); d=int(f[2])
    if lo+pad<=p<=hi-pad: ins.append(d)
    elif p<lo-pad or p>hi+pad: fl.append(d)
mi=statistics.median(ins) if ins else 0; mf=statistics.median(fl) if fl else 0
r=(mi/mf) if mf>0 else (float('inf') if mi>0 else float('nan'))
print('$v\t$c\t$lo\t$hi\t$sl\t$lib\t%.1f\t%.1f\t%.2f'%(mi,mf,r))" >> prophage_ratios.tsv
  done
done < embedded.tsv
echo PROPHAGE_DONE
