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
		* What are some errors you can spot in the data? Typos, missing data, inconsistencies?
		* What improvements could be made to the choices in naming?
		* What are some errors in the spreadsheet that would be difficult to spot? Is there anyway you can test this?
	* metadata back from facility: https://raw.githubusercontent.com/hbc/dc_2016_04/master/data/sequencing_results_metadata.txt
		* How are these samples organized?
		* If you wanted to relate file names to the sample names submitted above (e.g. wild type…) could you do so?
		* What do the _R1/_R2 extensions mean in the file names?
		* What does the ‘.gz’ extension on the filenames indicate?
		* What is the total file size - what challenges in downloading and sharing these data might exist?
	
## Cloud computing for genomics

* what is cloud computing?
	* on-demand, virtual computing (virtual machine, or cloud instance)
	* difference between remote and local
	* ask about previous experience
	* advantages: access, admin rights, preconfigured (BioLinux)
	* disadvantages: may cost $$$, may be difficult to obtain help
	* commercial: Amazon (using today), Google Cloud; open science: Atmosphere, JetStream
* logging in
	* Mac and Linux
		* terminal
		* ssh username@host
		* answer yes (first time)
		* enter password
	* Windows
		* putty
		* Host, dcuser, ps, port (22)
* data transfer
	* can download from URL using `wget URL`
	* can transfer using FileZilla (we'll do this later)
	* today the data and software we'll be using are already available on our custom cloud instance
	* quick download of some data files used today: http://www.datacarpentry.org/wrangling-genomics/variant_calling.tar.gz

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
	* make directories: ~/dc_workshop, containing docs, data, results
	* history

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

* link to module: http://www.datacarpentry.org/wrangling-genomics/
* required software (preinstalled): http://www.datacarpentry.org/wrangling-genomics/setup/

* instructor setup:
	* make sure data are back in hidden folder: `mv ~/.dc_sampledata_lite/untrimmed_fastq/ ~/dc_workshop/data/`
	* have two ssh sessions open to show things
	* have short list of notes available for students to reference?
	* draw directory structure on board
	* open lesson webpage to show data formats and other figures

* describe variant calling
	* data types
	* main process
	* all software we'll use is already available on the cloud instance
	* in many cases, there are multiple ways to accomplish the same task
		* some are more efficient (use less typing)
		* we may skip running steps on all files and copy processed data over to save time
		* we'll use absolute paths for this lesson to minimize confusion (but relative paths can work too!)

* Quality control
	* cd (home directory)
	* move data to project directory: `mv ~/.dc_sampledata_lite/untrimmed_fastq/ ~/dc_workshop/data/`
	* change directory to data folder: `cd ~/dc_workshop/data/untrimmed_fastq/`
	* software used in this section represents software downloaded by a user
		* executable files can be found wherever you decide to place them
	* **FastQC:** https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
	* run fastqc on untrimmed data: `~/FastQC/fastqc SRR097977.fastq`
	* can run on all untrimmed data with one command: `~/FastQC/fastqc *.fastq`
		* takes 4-5 minutes to run all (but don't need to run all)
		* in the meantime, discuss fasta/fastq format
	* directory for results: `mkdir ~/dc_workshop/results/fastqc_untrimmed_reads`
	* move results: `mv *.zip ~/dc_workshop/results/fastqc_untrimmed_reads/`, `mv *.html ~/dc_workshop/results/fastqc_untrimmed_reads/`
	* view results: change directory, try to unzip (`unzip *.zip`)
	* for loop: `for zip in *.zip; do unzip $zip; done`
	* save all results to file: `cat */summary.txt > ~/dc_workshop/docs/fastqc_summaries.txt`
	* **trimmomatic:** http://www.usadellab.org/cms/?page=trimmomatic
	* run with `java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar SE SRR098283.fastq \
SRR098283.fastq_trim.fastq SLIDINGWINDOW:4:20 MINLEN:20`
	* write in for loop: `for infile in *.fastq; do outfile=$infile\_trim.fastq; java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar SE $infile $outfile SLIDINGWINDOW:4:20 MINLEN:20; done`
		* takes 4-5 minutes to run all
		* in meantime, explain screen output, other options for trimming, can run FastQC on trimmed files to check
		* alternatively, copy already trimmed to directory: `mv ~/.dc_sampledata_lite/trimmed_fastq/ ~/dc_workshop/data/`
	* create new directory: `mkdir ~/dc_workshop/data/trimmmed_fastq`
	* move trimmed reads: `mv *_trim.fastq ~/dc_workshop/data/trimmed_fastq`

* Variant calling workflow
	* we'll work through one sequence sample, then talk about how to automate
	* software here is installed for computer, so path to executable does not need to be specified
	* start from `dc_workshop`: `cd ~/dc_workshop`
	* copy over reference genome data: `cp -r ~/.dc_sampledata_lite/ref_genome/ data/`
	* `mkdir results/sai results/sam results/bam results/bcf results/vcf`
	* **BWA:** http://bio-bwa.sourceforge.net/bwa.shtml
	* index reference: `bwa index data/ref_genome/ecoli_rel606.fasta`
	* align reads to reference: `bwa aln data/ref_genome/ecoli_rel606.fasta \
    data/trimmed_fastq/SRR097977.fastq_trim.fastq > results/sai/SRR097977.aligned.sai`
	* convert from SAI to SAM: `bwa samse data/ref_genome/ecoli_rel606.fasta \
    results/sai/SRR097977.aligned.sai \
    data/trimmed_fastq/SRR097977.fastq_trim.fastq > \
    results/sam/SRR097977.aligned.sam`
    * **samtools:** http://www.htslib.org/doc/samtools.html
    * convert SAM to BAM: `samtools view -S -b results/sam/SRR097977.aligned.sam > results/bam/SRR097977.aligned.bam`
    * sort BAM: `samtools sort results/bam/SRR097977.aligned.bam results/bam/SRR097977.aligned.sorted`
    * calculate read coverage: `samtools mpileup -g -f data/ref_genome/ecoli_rel606.fasta \
        results/bam/SRR097977.aligned.sorted.bam > results/bcf/SRR097977_raw.bcf`
    * **bcftools:** http://samtools.github.io/bcftools/bcftools.html
    * identify SNPs: `bcftools view -bvcg results/bcf/SRR097977_raw.bcf > results/bcf/SRR097977_variants.bcf`
    * filter SNPs: `bcftools view results/bcf/SRR097977_variants.bcf \ | /usr/share/samtools/vcfutils.pl varFilter - > results/vcf/SRR097977_final_variants.vcf`
    * overview vcf format: `less results/vcf/SRR097977_final_variants.vcf`
    * index BAM (required for visualization): `samtools index results/bam/SRR097977.aligned.sorted.bam`
    * transfer results files: 
    	* use FileZilla: enter Host, dcuser, ps, port (22) for Quickconnect
	* will take awhile to download (fairly large files), can demo IGV in meantime
	* orient to windows: select Desktop on left, navigate in cloud on right
	* data/ref_genome/ecoli_rel606.fasta, results/bam/SRR097977.aligned.sorted.bam, results/bam/SRR097977.aligned.sorted.bam.bai, results/vcf/SRR097977_final_variants.vcf
	* mention transferring data to a cloud instance
    
* Visualizing with IGV
	* Genomes/Load genomes from file: load reference genome
	* File/Load from file: load .bam (make sure .bai is in same location as .bam)
	* File/Load from file: load .vcf
	* Dark blue = heterozygous, Cyan = homozygous variant, Grey = reference
	* Filtered entries are transparent
	
* Automating with shell scripting
	* script available in dc_sample_data/variant_calling.tar.gz
	* uncompress and expand tarball: `tar -xvf variant_calling.tar.gz`

* continuing to work on your own
	* files (data and software) will be found in different places on different computers
	* the hardest part is keeping track of files, including upload/download
	* document your work thoroughly: easiest way is through scripts
