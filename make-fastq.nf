#!/usr/bin/env nextflow

nextflow.enable.dsl=2
BAM_GLOB = params.bam_glob


process sort_and_make_fastq {

    container 'docker.io/porchard/topmed_rnaseq_v2:v1'
    publishDir "${params.results}/fastq", mode: 'symlink'
    memory { 15.GB * task.attempt }
    time '72h'
    maxForks 100
    maxRetries 3
    errorStrategy {task.attempt <= maxRetries ? 'retry' : 'ignore'}
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-8]'
    cache 'lenient'
    tag "${tor}"

    input:
    tuple val(tor), path(bam)

    output:
    tuple val(tor), path("${tor}.1.fastq.gz"), path("${tor}.2.fastq.gz")

    """
    samtools sort -n -m 13G -o sorted.bam $bam
    java -Xmx${task.memory.getGiga() - 1}g -Xms${task.memory.getGiga() - 1}g -jar /opt/picard-tools/picard.jar SamToFastq I=sorted.bam FASTQ=${tor}.1.fastq.gz SECOND_END_FASTQ=${tor}.2.fastq.gz
    rm sorted.bam
    """

}

workflow {
    bams = Channel.fromPath(BAM_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // TOR, bam
    fastqs_sort = sort_and_make_fastq(bams)
}
