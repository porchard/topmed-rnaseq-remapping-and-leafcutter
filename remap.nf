#!/usr/bin/env nextflow

nextflow.enable.dsl=2
FASTQ_GLOB = params.fastq_glob
VCF_GLOB = params.vcf_glob
STAR_INDEX = params.star_index
TOR_NWD_MAPPINGS = params.tor_nwd_mappings

chroms = (1..22).collect({it -> "chr" + it}) + ['chrX']

// load TOR --> NWD mappings
TOR_TO_NWD = ['fake': '']
mapping_file = file(TOR_NWD_MAPPINGS)
for (line in mapping_file.readLines()) {
    x = line.replaceAll('\n', '').tokenize('\t')
    TOR_TO_NWD[x[0]] = x[1]
}


process filter_vcf {

    container 'library://porchard/default/general:20220107'
    memory '5 GB'
    clusterOptions='--partition=topmed-working,main --exclude=topmed,topmed2,topmed8'
    tag "${nwd}"
    maxForks 400
    publishDir "${params.results}/filtered-fastq"
    errorStrategy 'ignore'
    cache 'lenient'
    time '24h'

    input:
    tuple val(nwd), path("in.vcf.gz")

    output:
    tuple val(nwd), path("${nwd}.vcf.gz")

    """
    bcftools view -f 'PASS,.' -Oz -o ${nwd}.vcf.gz in.vcf.gz
    """

}

process remap_fastq {

    container 'docker.io/porchard/topmed_rnaseq_v2:v1'
    cpus 8
    memory '50 GB'
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed2,topmed8'
    tag "${tor}"
    maxForks 200
    publishDir "${params.results}/star"
    errorStrategy 'ignore'
    cache 'lenient'
    time '568h'

    input:
    tuple val(nwd), val(tor), path(fastq_1), path(fastq_2), path(vcf), path(star_index)

    output:
    tuple val(tor), path("${tor}.Aligned.sortedByCoord.out.bam"), path("${tor}.Aligned.sortedByCoord.out.bam.bai")

    """
    run_STAR.py --threads 8 --varVCFfile $vcf $star_index $fastq_1 $fastq_2 $tor
    """

}


process exon_exon_junction_counts {

    container '/net/snowwhite/home/porchard/singularity/leafcutter/2021-12-14/leafcutter_latest.sif'
    memory { 15.GB * task.attempt }
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed2,topmed[5-6],topmed8'
    tag "${tor}"
    maxForks 300
    maxRetries 2
    errorStrategy {task.attempt <= maxRetries ? 'retry' : 'ignore'}
    publishDir "${params.results}/exon-exon-junction-counts"
    cache 'lenient'
    
    input:
    tuple val(tor), path(bam), path(bam_index)

    output:
    tuple val(tor), path("${tor}.regtools_junc.txt.gz")

    """
    samtools view -h -q 255 $bam | grep -v "vW:i:[2-7]" | samtools view -b > filtered.bam
    samtools index filtered.bam
    regtools junctions extract -a 8 -m 50 -M 500000 -s 1 filtered.bam | gzip -c > ${tor}.regtools_junc.txt.gz
    rm filtered.bam filtered.bam.bai
    """

}

process filter_junction_file {

    clusterOptions='--partition=topmed-working --exclude=topmed,topmed[2-10]'
    tag "${tor}"
    maxForks 200
    maxRetries 2
    errorStrategy {task.attempt <= maxRetries ? 'retry' : 'ignore'}
    publishDir "${params.results}/exon-exon-junction-counts-filtered"
    cache 'lenient'

    input:
    tuple val(library), path('junctions.txt.gz')

    output:
    path("${library}.regtools_junc.txt.gz")

    """
    zcat junctions.txt.gz | grep -w -e ${chroms.join(' -e ')} | gzip -c > ${library}.regtools_junc.txt.gz
    """

}




workflow {
    first_fastq = Channel.fromPath(FASTQ_GLOB).filter({it -> it.getName().contains('.1.fastq.gz')}).map({it -> [it.getName().tokenize('.')[0], it]}) // TOR, fastq_1
    second_fastq = Channel.fromPath(FASTQ_GLOB).filter({it -> it.getName().contains('.2.fastq.gz')}).map({it -> [it.getName().tokenize('.')[0], it]}) // TOR, fastq_2
    fastq = first_fastq.combine(second_fastq, by: 0).filter({it -> TOR_TO_NWD.containsKey(it[0])}) // TOR, fastq1, fastq2
    vcfs = Channel.fromPath(VCF_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) | filter_vcf // NWD, vcf
    star_index = Channel.fromPath(STAR_INDEX)

    star_out = fastq.map({it -> [TOR_TO_NWD[it[0]]] + it}).combine(vcfs, by: 0).combine(star_index) | remap_fastq
    eejc = exon_exon_junction_counts(star_out) | filter_junction_file
}
