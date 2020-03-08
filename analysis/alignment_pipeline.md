# RNA-seq aligment


Aligment was perfomed using RSEM running with Bowtie2. RSEM is a software tool for quantifying transcript abundances from RNA-Seq data. RSEM mapped with Bowtie2 was used as it outperforms other aligment methods (Teng et al., 2016). 

    	module load rsem bowtie2 java-jdk/1.8.0_20


## Test quality control

	  for seqlib in $(cat minute_samples.txt)
	  do
	  fastqc raw_reads/fastq/${seqlib}*R1_001.fastq.gz -o raw_reads/fastqc_reports/
	  fastqc raw_reads/fastq/${seqlib}*R2_001.fastq.gz -o raw_reads/fastqc_reports/
	  done


The fastqc report shows a low amount of contaminants, adaptors, quality score, N content. It also shows that bases in the begining and end of the read have low quality, requiring trimming


## Trim adaptors and remove low quality reads (trimmonatic-0.38)

	  for seqlib in $(cat minute_samples.txt)
	  do
	  java -jar programs/Trimmomatic-0.38/Trimmomatic-0.38.jar PE -phred33 raw_reads/fastq/${seqlib}*R1_001.fastq.gz raw_reads/fastq/${seqlib}*R2_001.fastq.gz 	raw_reads/trimmed_reads/${seqlib}*R1_001.fastq.gz raw_reads/trimmed_reads/${seqlib}*R2_001.fastq.gz LEADING:3 TRAILING:3 MINLEN:36
	  done



## Unzip the reads

	  gunzip raw_reads/trimmed_reads/*.gz


## Download Drosophila genome


	wget ftp://ftp.ensembl.org/pub/release-94/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz

    	wget ftp://ftp.ensembl.org/pub/release-89/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.89.gtf.gz

## Alignment with RSEM


### Prepare reference for alignment

	  rsem-prepare-reference --gtf flyref/Drosophila_melanogaster.BDGP6.89.gtf --bowtie2 flyref/Drosophila_melanogaster.BDGP6.dna.toplevel.fa flyref/		Drosophila_melanogaster.BDGP5.89

### Alignment

	  for seqlib in $(cat minute_samples.txt)
	  do
	  rsem-calculate-expression --paired-end --bowtie2 --seed 123 raw_reads/trimmed_reads/${seqlib}*R1_001.fastq.gz raw_reads/trimmed_reads/${seqlib}	*R2_001.fastq.gz  flyref/RSEM_ref/Drosophila_melanogaster.BDGP5.70 rsem_results/${seqlib}
	  done 

Teng, M., Love, M.I., Davis, C.A., Djebali, S., Dobin, A., Graveley, B.R., Li, S., Mason, C.E., Olson, S., Pervouchine, D. and Sloan, C.A., 2016. A benchmark for RNA-seq quantification pipelines. Genome biology, 17(1), p.74.
