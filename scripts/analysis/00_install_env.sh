#!/bin/bash
set -euo pipefail
cd /NatureUsers/glitch/tools
[ -x micromamba ] || curl -sL -o micromamba https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
chmod +x micromamba
export MAMBA_ROOT_PREFIX=/NatureUsers/glitch/tools/mamba
./micromamba create -y -p /NatureUsers/glitch/envs/mvtools -c conda-forge -c bioconda \
  genomad phist minced megahit seqkit blast mmseqs2 python=3.11 pandas
echo INSTALL_DONE
