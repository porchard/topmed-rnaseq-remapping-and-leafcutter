#!/usr/bin/env nextflow

nextflow.enable.dsl=2

VCF_GLOB = params.vcf_glob
CHROMS = (1..22)
chroms = CHROMS.collect({it -> "chr" + it}) + ['chrX']


process list_samples_in_bcf {

    cache 'lenient'
    tag "${chrom}"
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed2,topmed[4-10]'

    input:
    tuple val(chrom), path(vcf)

    output:
    tuple val(chrom), path("sample-set.*")

    """
    bcftools view --header-only $vcf | grep CHROM | perl -pe 's/\\t/\\n/g' | awk 'NR>=10' > samples.txt
    split --lines 1000 samples.txt sample-set.
    """

}


process split_bcf {

    errorStrategy 'ignore'
    cache 'lenient'
    tag "${chrom}"
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed2,topmed[4-10]'

    input:
    tuple val(chrom), path(vcf), path(samples)

    output:
    tuple val(chrom), path("*.vcf.gz")

    """
    cat $samples | perl -pe 's/(.*)/\$1\\t-\\t\$1.${chrom}/' > samples-with-basenames.txt
    bcftools +split $vcf -Oz -o . -i'GT="alt"' -k FMT --samples-file samples-with-basenames.txt
    """

}


process merge_chroms {

    publishDir "${params.results}/vcf", mode: 'symlink'
    memory '2 GB'
    errorStrategy 'retry'
    cache 'lenient'
    tag "${subject_id}"
    maxForks 100
    clusterOptions='--partition=topmed-working --exclude=topmed,topmed2,topmed[4-10]'

    input:
    tuple val(subject_id), path(vcfs)

    output:
    tuple val(subject_id), path("${subject_id}.snps.vcf.gz")

    script:
    tabix_commands = vcfs.collect({x -> "tabix " + x}).join('\n')
    tmp = chroms.collect({x -> subject_id + "." + x + ".vcf.gz"}).join(' ')

    """
    $tabix_commands
    bcftools concat $tmp -Oz -o ${subject_id}.snps.vcf.gz
    rm ${subject_id}.chr*
    """

}



workflow {
    vcfs = Channel.fromPath(VCF_GLOB).map({it -> [it.getName().tokenize('.')[0], it]}) // chrom, vcf
    sample_lists = list_samples_in_bcf(vcfs).transpose() // chrom, sample_list
    split_bcf(vcfs.combine(sample_lists, by: 0)).transpose().map({it -> [it[1].getName().tokenize('.')[0], it[1]]}).groupTuple() | merge_chroms
}
