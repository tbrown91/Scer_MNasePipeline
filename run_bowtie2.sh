#!/bin/bash

bowtie2 -x ~/R_AN15mGEL/TOM_BROWN/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome -U BY_D_rep1.fastq -S ../sam_files/BY_D_rep1.sam &
bowtie2 -x ~/R_AN15mGEL/TOM_BROWN/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome -U BY_D_rep2.fastq -S ../sam_files/BY_D_rep2.sam &
bowtie2 -x ~/R_AN15mGEL/TOM_BROWN/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome -U BY_15mG_rep1.fastq -S ../sam_files/BY_15mG_rep1.sam &
bowtie2 -x ~/R_AN15mGEL/TOM_BROWN/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome -U BY_15mG_rep2.fastq -S ../sam_files/BY_15mG_rep2.sam &
bowtie2 -x ~/R_AN15mGEL/TOM_BROWN/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome -U BY_60mG_rep1.fastq -S ../sam_files/BY_60mG_rep1.sam &
bowtie2 -x ~/R_AN15mGEL/TOM_BROWN/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome -U BY_60mG_rep2.fastq -S ../sam_files/BY_60mG_rep2.sam &
