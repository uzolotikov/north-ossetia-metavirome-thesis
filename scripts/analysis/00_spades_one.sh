#!/bin/bash
s=$1; T=$2
E=/NatureUsers/nsinyavskiy/miniconda3/envs/metaspades-env/bin; C=/NatureUsers/nsinyavskiy/metaviroms_kgt/clean/pilot; O=/NatureUsers/glitch/mv_review/assembly_spades
$E/spades.py -1 $C/${s}_R1.fastq.gz -2 $C/${s}_R2.fastq.gz -o $O/$s -t $T -m 150 > $O/$s.log 2>&1 && echo "SPADES_ONE_DONE $s" || echo "SPADES_ONE_FAIL $s"
