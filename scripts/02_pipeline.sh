set -euo pipefail
cd "$(dirname "$0")/.."

THREADS=4
REF=ref/MN908947.3.fasta
BED=ref/nCoV-2019.primer.bed
NEXTCLADE_DATA=ref/nextclade_sars-cov-2
MIN_DEPTH=10            
MIN_FREQ=0.75           
MIN_QUAL=20             
mkdir -p results logs

# One-time preparation
[ -s "${REF}.bwt" ] || bwa index "$REF"
[ -s "${REF}.fai" ] || samtools faidx "$REF"
[ -d "$NEXTCLADE_DATA" ] || nextclade dataset get --name sars-cov-2 --output-dir "$NEXTCLADE_DATA"

process_sample() {
    local S=$1
    local D=results/$S
    mkdir -p "$D"
    if [ -s "$D/$S.consensus.fa" ]; then
        echo "$S: already processed, skipped"
        return
    fi

    echo "$S: 1/5 read QC"
    fastp -i "data/raw/${S}_1.fastq.gz" -I "data/raw/${S}_2.fastq.gz" \
          -o "$D/${S}_1.trim.fastq.gz" -O "$D/${S}_2.trim.fastq.gz" \
          --length_required 50 --thread "$THREADS" \
          --json "$D/fastp.json" --html "$D/fastp.html" > "$D/fastp.log" 2>&1

    echo "$S: 2/5 alignment"
    bwa mem -t "$THREADS" "$REF" "$D/${S}_1.trim.fastq.gz" "$D/${S}_2.trim.fastq.gz" 2> "$D/bwa.log" \
        | samtools sort -@ "$THREADS" -o "$D/$S.sorted.bam" - 2> /dev/null
    samtools index "$D/$S.sorted.bam"

    echo "$S: 3/5 primer trimming"
    ivar trim -i "$D/$S.sorted.bam" -b "$BED" -p "$D/$S.primertrim" -e -q "$MIN_QUAL" -m 30 \
        > "$D/ivar_trim.log" 2>&1                        # includes reads per primer (amplicon)
    samtools sort -@ "$THREADS" -o "$D/$S.primertrim.sorted.bam" "$D/$S.primertrim.bam" 2> /dev/null
    samtools index "$D/$S.primertrim.sorted.bam"

    echo "$S: 4/5 depth"
    samtools depth -a "$D/$S.primertrim.sorted.bam" > "$D/$S.depth.txt"

    echo "$S: 5/5 consensus"
    samtools mpileup -aa -A -d 0 -B -Q 0 "$D/$S.primertrim.sorted.bam" 2> /dev/null \
        | ivar consensus -p "$D/$S.consensus" -q "$MIN_QUAL" -t "$MIN_FREQ" -m "$MIN_DEPTH" -n N \
        > "$D/ivar_consensus.log" 2>&1
    sed -i "1s/.*/>$S/" "$D/$S.consensus.fa"             # header = sample name only

    rm -f "$D/${S}_1.trim.fastq.gz" "$D/${S}_2.trim.fastq.gz" "$D/$S.primertrim.bam"   # save space
}

while IFS=$'\t' read -r S lineage; do
    process_sample "$S" < /dev/null
done < data/samples.tsv

# Summary table: one row per sample
echo -e "sample\texpected_lineage\treads_raw\treads_after_qc\tpct_mapped\tmean_depth\tpct_depth_ge${MIN_DEPTH}\tn_bases" \
    > results/summary.tsv
while IFS=$'\t' read -r S lineage; do
    D=results/$S
    read -r raw qc < <(grep '"total_reads"' "$D/fastp.json" | head -2 | grep -o '[0-9]\+' | paste -sd ' ')
    mapped=$(samtools flagstat "$D/$S.sorted.bam" | grep "primary mapped" | grep -o '[0-9.]\+%' | tr -d '%')
    read -r mean pct < <(awk -v m="$MIN_DEPTH" '{ s += $3; if ($3 >= m) ok++ }
                          END { printf "%.0f %.2f\n", s/NR, 100*ok/NR }' "$D/$S.depth.txt")
    nb=$(grep -v ">" "$D/$S.consensus.fa" | tr -cd 'N' | wc -c)
    echo -e "$S\t$lineage\t$raw\t$qc\t$mapped\t$mean\t$pct\t$nb" >> results/summary.tsv
done < data/samples.tsv

# Lineages: all consensus genomes in one Nextclade run
cat results/*/*.consensus.fa > results/all_consensus.fasta
nextclade run --input-dataset "$NEXTCLADE_DATA" --output-tsv results/nextclade.tsv \
    results/all_consensus.fasta > logs/nextclade.log 2>&1

# Expected (benchmark, pangolin v3.1.3) vs Nextclade, columns found by name
awk -F'\t' -v OFS='\t' '
    NR == FNR { expected[$1] = $2; next }
    FNR == 1 {
        for (i = 1; i <= NF; i++) {
            if ($i == "seqName")          s = i
            if ($i == "clade")            c = i
            if ($i == "Nextclade_pango")  p = i
            if ($i == "qc.overallStatus") q = i
        }
        print "sample", "expected_lineage", "nextclade_pango", "nextclade_clade", "nextclade_qc"
        next
    }
    { print $s, expected[$s], (p ? $p : "NA"), $c, (q ? $q : "NA") }
' data/samples.tsv results/nextclade.tsv > results/lineages.tsv

echo
echo "=== Per-sample summary ==="
column -t -s $'\t' results/summary.tsv
echo
echo "=== Lineages: expected vs Nextclade ==="
column -t -s $'\t' results/lineages.tsv
