# Functional_impact_of_EMS-induced_mutations

The following repository containes the scripts used to identify EMS induced mutations, perform differential expression analysis, as well as the commands used to functionally annotate the mutations using SNPEff. The scripts are as follow:

The following BWA (Li and Durbin, 2010) command was used along with default parameters to align the fastq files of mutant lines to the D. pulicaria (Jackson et al. 2021) reference genomes:
bwa mem reference_genome.fa inputfile_1.fq.gz inputfile_2.fq.gz > output.sam

Samtools (Li et al. 2009) was used to remove reads that mapped to multiple locations and Picard tools's (http://broadinstitute.github.io/picard/) MarkDuplicates function was used to locate and tag PCR duplicates.
samtools view -h input.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b -o output.bam java -Xmx2g -XX:MaxMetaspaceSize=256m -jar picard.jar MarkDuplicates INPUT=input.bam OUTPUT=output.bam METRICS_FILE=file.metric

The mpileup and call functions of BCFtools (Li, 2011) along with default parameters were used to generate genotype likelihoods and genotype calls in a VCF file containing all EMS mutant lines derived from each natural Daphnia isolate. We added the following FORMAT and INFO tags to the VCF file: AD (allelic depth), DP (number of high-quality bases), ADF (allelic depth on forward strand) and ADR (allelic depth on reverse strand). For this study only biallelic single nucleotide polymorphisms (SNPs) that met the following filter parameters were used: Quality score (QUAL) >= 20, Sequencing depth (DP) >= 10, Distance >= 50 bp from an indel. Indels were not examined in this study and were also filtered out. An examle output file, output_file.vcf, can be viewed in this repository.
bcftools mpileup -Ou -f reference_genome.fa -a INFO/AD,FORMAT/AD,FORMAT/DP,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,FORMAT/SCR input_file_1.bam input_file_2.bam input_file_3.bam | bcftools call --threads xx -mO z -o output_file.vcf.gz

bcftools filter -S . -i 'QUAL>=20 & FORMAT/DP>=10 & FORMAT/DP<=60' input.vcf -Ov -o output.vcf ## Quality and depth filter

bcftools filter -g 50 input.vcf -Ov -o output.vcf ## Indel filter

After additional filtering with BCFtools, a custom python script, filter_mutations.py, was used to further filter mutations based on a consensus method, as well as forward and reverse read support. For each SNP site a consensus genotype call is established using a majority rule. For this study, 12 out of 16 samples needed to be in agreement to esablish a consensus call. If an EMS sample shows a genotype call different from the consensus, a tentative mutation is called. Since EMS has been shown to induce mutations randomly into the genome, it is highly unlikely that the same mutation will be observed amongst multiple EMS treatment lines. Each genotype call also has to be supported by at least 2 forward and 2 reverse reads in order to limit false positives due to sequencing errors. Both the number of samples needed for the consensus genotype and the read support can be manually adjusted as indicated in the script. The output files are a VCF file containing the final set of mutations, a text file containing the breakdown of the types of base substitutions, and a text file containing all of the mutations that failed the filtering criteria.
python3 filter_mutations.py

SNPEff was used to functionally annotate the final list of mutations. The -cancer option was used because it allows for the direct comparison between mutant genotypes and wildtype genotypes. The input file is the filtered mutations to be annotated and the output file contains the functionally annotated mutations. For more information on SNPEff please see https://pcingola.github.io/SnpEff/.
java -Xmx4g -jar snpEff.jar -v -cancer -cancerSamples untreated_treated_samples.txt ref_genome input.vcf > annotated_output.vcf

STAR aligner (Dobin et al. 2013) was used to map RNAseq reads to the D. pulicaria reference genome:
STAR --genomeDir raw_data \
--runThreadN 6 \
--readFilesIn forward_paired.fq.gz reverse_paired.fq.gz \
--readFilesCommand zcat \
--outFileNamePrefix prefix_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes Standard

FeatureCounts (Liao et al. 2014) was used to obtain the read counts for each gene:
featureCounts -p -T 8 -s 0 -a pulicaria.gtf -o sample.counts sample.bam

To identify differentially expressed genes, DEseq2 (Love et al. 2014) was used.
