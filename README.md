# bam2vcf - Variant identification based on bam file

Background:
A bioinformatics pipeline was designed to process targeted sequencing and WGS data, including several filters to eliminate non-MTBC reads and obtain reliable variant identification. The key analysis steps were read filtering based on the read alignment and variant filtering according to the base quality and allele depth.

Methods:
The Perl script bam2vcf detects variants using bam file with strict filtering conditions to minimize the impact of contamination and mixed infection (reads were filtered using the following criteria: read length<40, map quality<30, read quality<30, unpaired read, insert length>average insert length + 3SD, insert length<average insert length - 3SD, unproperly alignment read, soft clip read, hard clip read, mismatch>=5%, mapping identity<90%, or duplication read; variants were filtered by the following criteria: read depth <3, allelic depth <3 or minor allele frequency<1%).

Findings:
Using Illumina whole-exome data of HG002 and its high-confidence variant calls from GIAB as benchmark data, variants identified by our pipeline (called bam2vcf), gatk and varscan2 were compared. The F1 score, accuracy, and positive predictive value were higher in bam2vcf than in gatk and varscan2 when the alternative AF was ≥5% and alternative allele depth was ≥3X. 
