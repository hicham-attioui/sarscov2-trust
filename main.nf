params.samples        = "data/samples.tsv"
params.reads          = "data/raw"
params.ref            = "ref/MN908947.3.fasta"
params.primers        = "ref/nCoV-2019.primer.bed"
params.nextclade_data = "ref/nextclade_sars-cov-2"
params.outdir         = "results_nf"
params.min_depth      = 10       
params.min_freq       = 0.75     
params.min_qual       = 20       


process FASTP {
    tag "$sample"
    publishDir params.outdir, mode: 'copy', pattern: '*/fastp.*'

    input:
    tuple val(sample), val(scheme), path(r1), path(r2)

    output:
    tuple val(sample), val(scheme), path("${sample}_1.trim.fastq.gz"), path("${sample}_2.trim.fastq.gz"), emit: reads
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
    tuple path(ref), path("${ref}.*")      r

    script:
    """
    bwa index ${ref}
    """
}


process ALIGN {
    tag "$sample"
    publishDir params.outdir, mode: 'copy', pattern: '*/flagstat.txt'

    input:
    tuple val(sample), val(scheme), path(r1), path(r2)
    tuple path(ref), path(index)           

    output:
    tuple val(sample), val(scheme), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam
    path "${sample}/flagstat.txt"

    script:
    """
    mkdir -p ${sample}
    bwa mem -t ${task.cpus} ${ref} ${r1} ${r2} 2> bwa.log \\
        | samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam -
    samtools index ${sample}.sorted.bam
    samtools flagstat ${sample}.sorted.bam > ${sample}/flagstat.txt
    """
}


process TRIM_PRIMERS {
    tag "$sample"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample), val(scheme), path(bam), path(bai)
    path primers

    output:
    tuple val(sample), path("${sample}/${sample}.primertrim.sorted.bam"),
                       path("${sample}/${sample}.primertrim.sorted.bam.bai"), emit: bam
    path "${sample}/ivar_trim.log"

    script:
    if (scheme == 'ARTIC_V3')
        """
        mkdir -p ${sample}
        ivar trim -i ${bam} -b ${primers} -p ${sample}.primertrim -e -q ${params.min_qual} -m 30 \\
            > ${sample}/ivar_trim.log 2>&1
        samtools sort -@ ${task.cpus} -o ${sample}/${sample}.primertrim.sorted.bam ${sample}.primertrim.bam
        samtools index ${sample}/${sample}.primertrim.sorted.bam
        """
    else
        """
        mkdir -p ${sample}
        echo "Primer scheme ${scheme}: ARTIC V3 primer trimming skipped" > ${sample}/ivar_trim.log
        cp ${bam} ${sample}/${sample}.primertrim.sorted.bam
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
    samtools depth -aa ${bam} > ${sample}/${sample}.depth.txt
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
    path consensus                        
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
    reads = channel.fromPath(params.samples)
        .splitCsv(sep: '\t')
        .map { row ->
            def scheme = row.size() >= 3 ? row[1] : 'ARTIC_V3'
            tuple(row[0], scheme,
                  file("${params.reads}/${row[0]}_1.fastq.gz", checkIfExists: true),
                  file("${params.reads}/${row[0]}_2.fastq.gz", checkIfExists: true)) }

    BWA_INDEX(file(params.ref, checkIfExists: true))
    FASTP(reads)
    ALIGN(FASTP.out.reads, BWA_INDEX.out)
    TRIM_PRIMERS(ALIGN.out.bam, file(params.primers, checkIfExists: true))
    DEPTH(TRIM_PRIMERS.out.bam)
    CONSENSUS(TRIM_PRIMERS.out.bam)
    NEXTCLADE(CONSENSUS.out.collect(sort: true), file(params.nextclade_data, checkIfExists: true))
}
