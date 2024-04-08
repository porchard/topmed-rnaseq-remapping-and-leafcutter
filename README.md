# WASP remapping of TOPMed RNA-seq data

## Dependencies

* Singularity (v. >=3)
* NextFlow (v. >= 21.04.0)


## Setup

1. Pull a singularity container with some dependencies: `make singularity`
2. Fetch data: `make data`
3. Place BAM files (named: `{TOR}.bam` into `data/bam`)
4. Place BCF files (named: `{chrom}.bcf`, containing only the samples of interest) into `data/bcf`
5. Place a file named 'TOR-to-NWD.txt' in `data/tor-to-nwd`. Should be a two column TSV file (no header) with RNA-seq ID (TOR) in the first column and WGS ID (NWD) in the second. These are the samples that will be remapped. The corresponding RNA-seq BAM file (`{TOR}.bam`) should be in `data/bam`, and there should be a sample with that NWD ID in the BCF files.


## Running

1. Convert BAM files to fastq: `make fastq`
2. Extract per-sample VCF files: `make fetch-genotypes`
3. Remap: `make remap`