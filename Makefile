ROOT=/net/topmed11/working/porchard/rnaseq-wasp-remap
DATA=$(ROOT)/data
WORK=$(ROOT)/work
BIN=$(ROOT)/bin

ANALYSIS=$(WORK)/$@

.PHONY: all

define NL


endef

##### DATA #####
data: bam bcf tor-to-nwd star-index

star-index:
	mkdir -p $(DATA)/$@ && wget https://personal.broadinstitute.org/francois/topmed/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v30_oh100.tar.gz && tar -xvzf STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v30_oh100.tar.gz

bam:
	mkdir -p $(DATA)/$@
	ln -s /net/topmed10/working/porchard/rnaseq/data/tor-bam/* $(DATA)/$@/
	ln -s /net/topmed10/working/porchard/rnaseq/data/tor-bam-mesa/* $(DATA)/$@/

bcf:
	mkdir -p $(DATA)/$@
	ln -s /net/topmed10/working/porchard/rnaseq/work/subset-topmed-bcf/freeze-alpha/results/bcfs-by-chrom/chr* $(DATA)/$@/

tor-to-nwd:
	mkdir -p $(DATA)/$@
	cp /net/topmed10/working/porchard/rnaseq/work/remap-for-sqtl-pass-only/data/tor-to-nwd.txt $(DATA)/$@/

##### ANALYSES #####
fastq:
	cd $(ANALYSIS) && nohup nextflow run -resume -queue-size 100 --bam_glob '$(DATA)/bam/*.bam' --results $(ANALYSIS)/results $(ROOT)/make-fastq.nf &

# Note: phasing / no phasing shouldn't matter: https://github.com/alexdobin/STAR/issues/403
fetch-genotypes:
	mkdir -p $(ANALYSIS)/data/vcf
	cd $(ANALYSIS) && nohup nextflow run -resume --vcf_glob '$(DATA)/bcf/*' --results $(ANALYSIS)/results $(ROOT)/$@.nf &

remap:
	mkdir -p $(ANALYSIS)/data/fastq
	cd $(ANALYSIS) && nohup nextflow run -resume -qs 1000 --results $(ANALYSIS)/results --fastq_glob '$(WORK)/fastq/results/fastq/TOR*.fastq.gz' --vcf_glob '$(WORK)/fetch-genotypes/results/vcf/*' --star_index $(DATA)/star-index/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v30_oh100 --tor_nwd_mappings $(DATA)/tor-to-nwd/TOR-to-NWD.txt $(ROOT)/remap.nf &