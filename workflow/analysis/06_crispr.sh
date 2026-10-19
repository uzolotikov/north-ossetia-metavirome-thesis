#!/bin/bash
# Спейсеры CRISPR (MinCED, >=3 повторов) во всех контигах сборок фагомов, клеточных сборок Тиб и сборки
# Бирагзанга, поиск их в каталоге vOTU: совпадение по всей длине спейсера, не более одного несовпадения.
source "$(dirname "$0")/../common.sh"
set +e
X=$ENVS/mvtools/bin; export PATH=$X:$PATH
W=$ANALYSIS
mkdir -p $W/hosts/crispr && cd $W/hosts/crispr
# все контиги, а не только >=5 кб: массивы CRISPR часто лежат на коротких контигах
for a in $MV; do ln -sf $OUT/assembly/no_NC/$a/contigs.fasta $a.fna; done
for a in $TIB_CELLULAR; do ln -sf $OUT/assembly/cellular/$a/contigs.fasta cell-$a.fna; done
ln -sf $BIRAGZANG Biragzang.fna
for f in *.fna; do b=${f%.fna}; $X/minced -minNR 3 -spacers $f $b.crisprs $b.gff > /dev/null 2>&1; [ -f ${b}_spacers.fa ] && sed "s/^>/>$b|/" ${b}_spacers.fa; done > all_spacers.fa
echo "arrays per assembly:"; for f in *.gff; do echo "$f $(grep -vc '^#' $f)"; done
echo "spacers: $(grep -c '>' all_spacers.fa)"
$X/blastn -task blastn-short -query all_spacers.fa -subject $OUT/votu/no_NC/catalog/votu_no_NC_95_85/representatives_ge5kb.fna -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen sstart send evalue" -evalue 1 -word_size 7 -dust no > spacers_vs_votu.raw.tsv
awk -F"\t" '$4==$7 && ($5+$6)<=1' spacers_vs_votu.raw.tsv > spacers_vs_votu.hits.tsv
echo "spacer hits (full length, <=1 mismatch): $(wc -l < spacers_vs_votu.hits.tsv)"
echo CRISPR_DONE
