#!/usr/bin/env bash
# Download the 24 failed-QC benchmark samples listed in data/samples_failedqc.tsv.
# Same logic as 01_download.sh: --split-3, unpaired reads removed, pairs checked,
# re-paired with seqkit when the two files differ. Safe to rerun (finished samples skipped).
set -uo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/raw logs
i=0; total=$(wc -l < data/samples_failedqc.tsv)
while IFS=$'\t' read -r run scheme desc; do
  i=$((i + 1))
  if [ -s data/raw/${run}_2.fastq.gz ]; then echo "[$i/$total] $run: already there"; continue; fi
  rm -f data/raw/${run}*
  echo "[$i/$total] $run: downloading ($desc)"
  fasterq-dump --split-3 --threads 4 --temp data/raw --outdir data/raw "$run" < /dev/null >> logs/download_failedqc.log 2>&1
  rm -f data/raw/${run}.fastq
  if [ ! -s data/raw/${run}_2.fastq ]; then echo "  WARNING: $run has no read 2 (not paired)"; continue; fi
  gzip -f data/raw/${run}_1.fastq data/raw/${run}_2.fastq
  n1=$(( $(zcat data/raw/${run}_1.fastq.gz | wc -l) / 4 )); n2=$(( $(zcat data/raw/${run}_2.fastq.gz | wc -l) / 4 ))
  if [ "$n1" -ne "$n2" ]; then
    seqkit pair -j 4 -1 data/raw/${run}_1.fastq.gz -2 data/raw/${run}_2.fastq.gz -O data/raw/paired 2> /dev/null
    mv data/raw/paired/${run}_?.fastq.gz data/raw/
    echo "  re-paired ($n1 / $n2)"
  fi
done < data/samples_failedqc.tsv
echo "Done."
