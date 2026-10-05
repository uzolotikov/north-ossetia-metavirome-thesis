#!/bin/bash
set -uo pipefail
export PATH=/NatureUsers/glitch/envs/mvtools/bin:$PATH
W=/NatureUsers/glitch/mv_review; X=/NatureUsers/glitch/envs/mvtools/bin; K=/NatureUsers/nsinyavskiy/metaviroms_kgt
mkdir -p $W/hosts/crispr && cd $W/hosts/crispr
# all contigs (not only >=5kb) of every assembly: CRISPR arrays often sit on short contigs
for a in MV-4122_S12 MV-4125_S18 MV-4128_S29; do ln -sf $K/assembly/no_NC/$a/contigs.fasta $a.fna; done
for a in 4122-0-22_S11 4122-0-44_S10 4125-0-22_S17 4125-0-44_S16; do ln -sf $K/assembly/cellular/$a/contigs.fasta cell-$a.fna; done
ln -sf /NatureUsers/oscypek/data/alania-meta/Biragzang.final.fna Biragzang.fna
for f in *.fna; do b=${f%.fna}; $X/minced -minNR 3 -spacers $f $b.crisprs $b.gff > /dev/null 2>&1; [ -f ${b}_spacers.fa ] && sed "s/^>/>$b|/" ${b}_spacers.fa; done > all_spacers.fa
echo "arrays per assembly:"; for f in *.gff; do echo "$f $(grep -vc '^#' $f)"; done
echo "spacers: $(grep -c '>' all_spacers.fa)"
$X/blastn -task blastn-short -query all_spacers.fa -subject $K/votu/no_NC/catalog/votu_no_NC_95_85/representatives_ge5kb.fna -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen sstart send evalue" -evalue 1 -word_size 7 -dust no > spacers_vs_votu.raw.tsv
awk -F"\t" '$4==$7 && ($5+$6)<=1' spacers_vs_votu.raw.tsv > spacers_vs_votu.hits.tsv
echo "spacer hits (full length, <=1 mismatch): $(wc -l < spacers_vs_votu.hits.tsv)"
echo CRISPR_DONE
