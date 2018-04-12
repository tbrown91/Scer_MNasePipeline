#!/bin/bash

samtools view -bS -F 3844 BY_D_rep1.sam > ../bam_files/BY_D_rep1.bam &
samtools view -bS -F 3844 BY_D_rep2.sam > ../bam_files/BY_D_rep2.bam &
samtools view -bS -F 3844 BY_15mG_rep1.sam > ../bam_files/BY_15mG_rep1.bam &
samtools view -bS -F 3844 BY_15mG_rep2.sam > ../bam_files/BY_15mG_rep2.bam &
samtools view -bS -F 3844 BY_60mG_rep1.sam > ../bam_files/BY_60mG_rep1.bam &
samtools view -bS -F 3844 BY_60mG_rep2.sam > ../bam_files/BY_60mG_rep2.bam &
