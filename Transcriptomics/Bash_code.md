# 1. Reference Data Download and Preparation

We download and decompress *Arabidopsis thaliana* (Arabidopsis)
reference genome and annotation files.

    wget -O A_thaliana.fa.gz https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-61/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz 

    wget -O A_thaliana.gtf.gz https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-61/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.61.gtf.gz -O A_thaliana.gtf

    gunzip A_thaliana.fa.gz
    gunzip A_thaliana.gtf.gz

# 2. Genome Index Generation

We build a STAR genome index using the Arabidopsis reference genome and
GTF file.

    STAR --runMode genomeGenerate \
         --genomeDir /home/omicas/transcriptomics/genome/index \
         --genomeFastaFiles A_thaliana.fa \
         --sjdbGTFfile A_thaliana.gtf \
         --genomeSAindexNbases 10

# 3. Paired-End Sample Processing Script

We define a bash script called “new\_sampr\_paired\_end.sh” for
paired-end samples processing.

    #! /bin/bash

    # Made by Aiub Mohamed Barara to process A.thaliana samples from RNA-seq 

    ## Script generalization
    FOLDER=$1
    NAME=$2

    ## Access to samples folder
    cd $FOLDER

    ## Sample quality control and read mapping to reference genome
    fastqc ${NAME}.R1.fastq.gz
    fastqc ${NAME}.R2.fastq.gz

    ## Read mapping
    STAR --genomeDir ../../genome/index/\
         --readFilesIn ${NAME}.R1.fastq.gz ${NAME}.R2.fastq.gz\
         --readFilesCommand 'gunzip -c' --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif\
         --outFilterIntronMotifs RemoveNoncanonical --alignIntronMax 300000\
         --outFileNamePrefix $NAME

    ## Transcript assembly and Quantification
    stringtie -G ../../genome/A_thaliana.gtf -o $NAME.gtf\
           ${NAME}Aligned.sortedByCoord.out.bam

    stringtie -e -G ../../genome/A_thaliana.gtf\
          -o $NAME.gtf -A $NAME.tsv\
          ${NAME}Aligned.sortedByCoord.out.bam

# 4. Single-End Sample Processing Script

We define a bash script called “new\_sampr\_single.sh” for single-end
samples processing.

    #! /bin/bash

    # Made by Aiub Mohamed Barara to process A.thaliana samples from RNA-seq 

    ## Script generalization
    FOLDER=$1
    NAME=$2

    ## Access to samples folder
    cd $FOLDER

    ## Sample quality control and read mapping to reference genome
    fastqc ${NAME}.fastq.gz

    ## Read mapping
    STAR --genomeDir ../../genome/index/\
         --readFilesIn ${NAME}.fastq.gz\
         --readFilesCommand 'gunzip -c' --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif\
         --outFilterIntronMotifs RemoveNoncanonical --alignIntronMax 300000\
         --outFileNamePrefix $NAME

    ## Transcript assembly and Quantification
    stringtie -G ../../genome/A_thaliana.gtf -o $NAME.gtf\
           ${NAME}Aligned.sortedByCoord.out.bam

    stringtie -e -G ../../genome/A_thaliana.gtf\
          -o $NAME.gtf -A $NAME.tsv\
          ${NAME}Aligned.sortedByCoord.out.bam

# 5. Job Submission via SLURM

We submit paired-end (such as “ntrc\_adult\_1.fq.gz”) and single-end
(such as “wt\_seed\_3.fq.gz”) sample processing jobs to the cluster.
Each command runs one of the two previously defined scripts: one for
paired-end (new\_sampr\_paired\_end.sh) and another for single-end
(new\_sampr\_single.sh) data.

    # Example for a paired-end sample
    sbatch --job-name=ntrc_adult_1 --output=ntrc_adult_1 new_sampr_paired_end.sh /home/omicas/transcriptomics/samples/NTRC_adult_1/ ntrc_adult_1

    # Example for a single-end sample
    sbatch --job-name=wt_seed_3 --output=wt_seed_3 new_sampr_single.sh /home/omicas/transcriptomics/samples/WT_seed_3/ wt_seed_3

# 6. Gene-Level Count Matrix Generation

We generate a count matrix (.csv file) for discrete expression analysis
using the prepDE.py script.

    cd /home/omicas/transcriptomics/
    prepDE.py -i samples/
