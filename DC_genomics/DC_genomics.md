## Required software

* Java
* Filezilla
* Integrative Genomics Viewer (IGV)
* Putty (for Windows only)

## Introduction to the data

* Introduction
	* dataset: http://www.datacarpentry.org/organization-genomics/01-introduction/
	* SRA: http://www.ncbi.nlm.nih.gov/sra?term=SRA026813
	* Click on first entry (ZDB30), All runs, RunInfo Table
	* downloads as SraRunTable.txt

* Planning for NGS projects
	* sample submission sheet: https://raw.githubusercontent.com/hbc/dc_2016_04/master/data/sample_submission.txt
	* metadata back from facility: https://raw.githubusercontent.com/hbc/dc_2016_04/master/data/sequencing_results_metadata.txt
	
## Cloud computing for genomics

* what is cloud computing?
* logging in
* file transfer: Filezilla and wget
* download data files
	* http://www.datacarpentry.org/wrangling-genomics/variant_calling.tar.gz

## Introduction to the command line/shell

* Introduction and the file system
	* why use the shell?
	* how to access the shell	
	* files and directories (slashes)
	* home and root directories
	* full and relative paths
	* cd
	* ls
	* man and --help
	* tab completion

* Project organization
	* mkdir
	* history
	* tail
	* | (pipe)
	* nano

* Working with files
	* cat
	* head
	* tail
	* cp
	* mv
	* rm

* Redirection
	* grep (use -B1 -A2 for line before and two lines after)
	* `>` and `>>`
	* cut
	* sort

* Writing scripts
	* creating .sh 
	* chmod

## Data wrangling and processing

* describe variant calling
	* data types
	* main process

* Quality control
	* cd (home directory)
	* mkdir dc_workshop
	* mkdir dc_workshop/data dc_workshop/docs dc_workshop/results
	* move data to project directory
	* change directory to data folder
	* run fastqc on untrimmed data
	* view results: change directory, try to unzip
	* for loop: for zip in *.zip; do unzip $zip; done
	* cat */summary.txt > ~/dc_workshop/docs/fastqc_summaries.txt
	* trimmomatic: run with java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar inputfile outputfile OPTION:VALUE...
	* write in for loop: for infile in *.fastq; do outfile=$infile\_trim.fastq; java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar SE $infile $outfile SLIDINGWINDOW:4:20 MINLEN:20; done

* Variant calling workflow
	* mkdir  results/sai results/sam results/bam results/bcf results/vcf
	* bwa index data/ref_genome/ecoli_rel606.fasta
	* bwa aln path/to/ref_genome.fasta path/to/fastq > SAIfile
	* bwa samse data/ref_genome/ecoli_rel606.fasta \
    results/sai/SRR097977.aligned.sai \
    data/trimmed_fastq/SRR097977.fastq_trim.fastq > \
    results/sam/SRR097977.aligned.sam
    * samtools view -S -b results/sam/SRR097977.aligned.sam > results/bam/SRR097977.aligned.bam
    * samtools sort results/bam/SRR097977.aligned.bam results/bam/SRR097977.aligned.sorted
    * samtools mpileup -g -f data/ref_genome/ecoli_rel606.fasta \
        results/bam/SRR097977.aligned.sorted.bam > results/bcf/SRR097977_raw.bcf
    * bcftools view -bvcg results/bcf/SRR097977_raw.bcf > results/bcf/SRR097977_variants.bcf
    * bcftools view results/bcf/SRR097977_variants.bcf \ | /usr/share/samtools/vcfutils.pl varFilter - > results/vcf/SRR097977_final_variants.vcf
    * overview vcf format
    * samtools index results/bam/SRR097977.aligned.sorted.bam
    * transfer results files: results/bam/SRR097977.aligned.sorted.bam, results/bam/SRR097977.aligned.sorted.bam.bai, data/ref_genome/ecoli_rel606.fasta, results/vcf/SRR097977_final_variants.vcf
    
* Visualizing with IGV
	* Genomes/Load genomes from file: load reference genome
	* File/Load from file: load .bam (make sure .bai is in same location as .bam)
	* File/Load from file: load .vcf
	* Dark blue = heterozygous, Cyan = homozygous variant, Grey = reference
	* Filtered entries are transparent
	
* Automating with shell scripting
