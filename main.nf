#!/usr/bin/env nextflow
/*
 * SARS-CoV-2 amplicon pipeline (ARTIC V3, Illumina paired-end).
 * Nextflow version of scripts/02_pipeline.sh: same tools, same commands, same parameters,
 * so the consensus genomes must be identical to the Bash version.
 *
 *   FASTP        read QC
 *   BWA_INDEX    index the reference (once)
 *   ALIGN        bwa mem + sort + index
 *   TRIM_PRIMERS ivar trim + sort + index (log keeps reads per primer)
 *   DEPTH        samtools depth -a
 *   CONSENSUS    samtools mpileup | ivar consensus
 *   NEXTCLADE    clade and Pango lineage of all genomes
 *
 * Per-sample results are written by each task into a folder named after the sample,
 * declared as outputs, and copied by publishDir to results_nf/<sample>/.
 * (publishDir itself is read before any sample exists, so it cannot use the sample name.)
 *
 * Run from the repository root, in the sarscov2 environment:
 *   nextflow run main.nf                 # laptop (local executor)
 *   nextflow run main.nf -resume         # reuse finished tasks
 *   nextflow run main.nf -profile slurm  # HPC cluster (see nextflow.config)
 */

params.samples        = "data/samples.tsv"
params.reads          = "data/raw"
params.ref            = "ref/MN908947.3.fasta"
params.primers        = "ref/nCoV-2019.primer.bed"
params.nextclade_data = "ref/nextclade_sars-cov-2"
params.outdir         = "results_nf"
params.min_depth      = 10       // minimum reads to call a base
params.min_freq       = 0.75     // minimum share of reads agreeing on the base
params.min_qual       = 20       // minimum base quality


process FASTP {
    tag "$sample"
    publishDir params.outdir, mode: 'copy', pattern: '*/fastp.*'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_1.trim.fastq.gz"), path("${sample}_2.trim.fastq.gz"), emit: reads
    path "${sample}/fastp.*"

    script:
    """
    mkdir -p ${sample}
    fastp -i ${r1} -I ${r2} \\
          -o ${sample}_1.trim.fastq.gz -O ${sample}_2.trim.fastq.gz \\
          --length_required 50 --thread ${task.cpus} \\
          --json ${sample}/fastp.json --html ${sample}/fastp.html
    """
}


process BWA_INDEX {
    input:
    path ref

    output:
    tuple path(ref), path("${ref}.*")      // the reference and its index files, kept together

    script:
    """
    bwa index ${ref}
    """
}


process ALIGN {
    tag "$sample"

    input:
    tuple val(sample), path(r1), path(r2)
    tuple path(ref), path(index)           // value channel: reused by every sample

    output:
    tuple val(sample), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai")

    script:
    """
    bwa mem -t ${task.cpus} ${ref} ${r1} ${r2} 2> bwa.log \\
        | samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam -
    samtools index ${sample}.sorted.bam
    """
}


process TRIM_PRIMERS {
    tag "$sample"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)
    path primers

    output:
    tuple val(sample), path("${sample}/${sample}.primertrim.sorted.bam"),
                       path("${sample}/${sample}.primertrim.sorted.bam.bai"), emit: bam
    path "${sample}/ivar_trim.log"

    script:
    """
    mkdir -p ${sample}
    ivar trim -i ${bam} -b ${primers} -p ${sample}.primertrim -e -q ${params.min_qual} -m 30 \\
        > ${sample}/ivar_trim.log 2>&1
    samtools sort -@ ${task.cpus} -o ${sample}/${sample}.primertrim.sorted.bam ${sample}.primertrim.bam
    samtools index ${sample}/${sample}.primertrim.sorted.bam
    """
}


process DEPTH {
    tag "$sample"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path "${sample}/${sample}.depth.txt"

    script:
    """
    mkdir -p ${sample}
    samtools depth -a ${bam} > ${sample}/${sample}.depth.txt
    """
}


process CONSENSUS {
    tag "$sample"
    publishDir params.outdir, mode: 'copy', pattern: '*/*.consensus.fa'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path "${sample}/${sample}.consensus.fa"

    script:
    """
    mkdir -p ${sample}
    samtools mpileup -aa -A -d 0 -B -Q 0 ${bam} 2> /dev/null \\
        | ivar consensus -p ${sample}/${sample}.consensus -q ${params.min_qual} -t ${params.min_freq} \\
                         -m ${params.min_depth} -n N
    sed -i "1s/.*/>${sample}/" ${sample}/${sample}.consensus.fa
    """
}


process NEXTCLADE {
    publishDir params.outdir, mode: 'copy'

    input:
    path consensus                         // all genomes at once (collected)
    path dataset

    output:
    path "all_consensus.fasta"
    path "nextclade.tsv"

    script:
    """
    cat ${consensus} > all_consensus.fasta
    nextclade run --input-dataset ${dataset} --output-tsv nextclade.tsv all_consensus.fasta
    """
}


workflow {
    // one item per sample: (sample, read 1, read 2), from the same sample list as the Bash version
    reads = channel.fromPath(params.samples)
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0],
                            file("${params.reads}/${row[0]}_1.fastq.gz", checkIfExists: true),
                            file("${params.reads}/${row[0]}_2.fastq.gz", checkIfExists: true)) }

    BWA_INDEX(file(params.ref, checkIfExists: true))
    FASTP(reads)
    ALIGN(FASTP.out.reads, BWA_INDEX.out)
    TRIM_PRIMERS(ALIGN.out, file(params.primers, checkIfExists: true))
    DEPTH(TRIM_PRIMERS.out.bam)
    CONSENSUS(TRIM_PRIMERS.out.bam)
    NEXTCLADE(CONSENSUS.out.collect(sort: true), file(params.nextclade_data, checkIfExists: true))
}
