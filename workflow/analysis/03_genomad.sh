#!/bin/bash
# geNomad 1.12.0 (база 1.9) по каталогу vOTU и по общему референсу step6. База скачивается в $DB/genomad_db.
source "$(dirname "$0")/../common.sh"
set +e
X=$ENVS/mvtools/bin; export PATH=$X:$PATH
W=$ANALYSIS
cd $DB
if [ ! -d genomad_db ]; then
  curl -sSL -o genomad_db_v1.9.tar.gz "https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz?download=1" && tar xzf genomad_db_v1.9.tar.gz && rm -f genomad_db_v1.9.tar.gz && echo GENOMAD_DB_OK || { echo GENOMAD_DB_FAIL; exit 1; }
fi
genomad end-to-end --cleanup --splits 4 -t $THREADS $OUT/votu/no_NC/catalog/votu_no_NC_95_85/representatives_ge5kb.fna $W/genomad/votu97 $DB/genomad_db > $W/logs/genomad_votu97.log 2>&1 && echo GENOMAD_VOTU_DONE || echo GENOMAD_VOTU_FAIL
genomad end-to-end --cleanup --splits 8 -t $THREADS $W/step6/ref_derep.fna $W/genomad $DB/genomad_db > $W/logs/genomad_ref.log 2>&1 && echo GENOMAD_REF_DONE || echo GENOMAD_REF_FAIL
