set -euo pipefail                 
cd "$(dirname "$0")/.."            

THREADS=4
mkdir -p data/raw ref logs

# 1. Reference genome
if [ ! -s ref/MN908947.3.fasta ]; then
    echo "Downloading reference genome..."
    wget -q -O ref/MN908947.3.fasta \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MN908947.3&rettype=fasta"
fi

# 2. ARTIC V3 primer scheme (BED file: one line per primer, with its start and end on the reference)
if [ ! -s ref/nCoV-2019.primer.bed ]; then
    echo "Downloading ARTIC V3 primer scheme..."
    wget -q -O ref/nCoV-2019.primer.bed \
        "https://github.com/artic-network/artic-ncov2019/raw/master/primer_schemes/nCoV-2019/V3/nCoV-2019.primer.bed"
fi

# 3. Benchmark table -> data/samples.tsv (run accession and expected lineage)
if [ ! -d data/benchmark ]; then
    git clone -q https://github.com/CDCgov/datasets-sars-cov-2.git data/benchmark
fi
TSV=data/benchmark/datasets/sars-cov-2-voivoc.tsv

awk -F'\t' '
    /^GISAID_assession/ {
        for (i = 1; i <= NF; i++) {
            if ($i == "SRArun_acc") run_col = i
            if ($i ~ /^Lineage/)    lin_col = i
        }
        next
    }
    run_col && $run_col ~ /^[SE]RR/ { print $run_col "\t" $lin_col }
' "$TSV" > data/samples.tsv

echo "Samples to download: $(wc -l < data/samples.tsv)"

# 4. Download each sample (paired-end reads), then compress
while IFS=$'\t' read -r run lineage; do
    if [ -s "data/raw/${run}_1.fastq.gz" ] && [ -s "data/raw/${run}_2.fastq.gz" ]; then
        echo "$run ($lineage): already downloaded, skipped"
        continue
    fi
    echo "$run ($lineage): downloading..."
    fasterq-dump --split-3 --threads "$THREADS" --temp data/raw --outdir data/raw "$run" \
        < /dev/null >> logs/download.log 2>&1
    rm -f "data/raw/${run}.fastq"                        
    gzip -f "data/raw/${run}_1.fastq" "data/raw/${run}_2.fastq"
done < data/samples.tsv

# Check: both files of a pair must contain the same number of reads
echo
echo "Pair check (reads in _1 and _2 must be equal):"
while IFS=$'\t' read -r run lineage; do
    n1=$(( $(zcat "data/raw/${run}_1.fastq.gz" | wc -l) / 4 ))
    n2=$(( $(zcat "data/raw/${run}_2.fastq.gz" | wc -l) / 4 ))
    if [ "$n1" -eq "$n2" ]; then status="OK"; else status="MISMATCH"; fi
    echo "  $run: $n1 / $n2  $status"
done < data/samples.tsv

echo
echo "Done. Files in data/raw:"
ls -lh data/raw
