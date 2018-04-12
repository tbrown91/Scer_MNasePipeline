#!/bin/bash

samtools sort BY_D_rep1.bam BY_D_rep1_sorted &
samtools sort BY_D_rep2.bam BY_D_rep2_sorted &
samtools sort BY_15mG_rep1.bam BY_15mG_rep1_sorted &
samtools sort BY_15mG_rep2.bam BY_15mG_rep2_sorted &
samtools sort BY_60mG_rep1.bam BY_60mG_rep1_sorted &
samtools sort BY_60mG_rep2.bam BY_60mG_rep2_sorted &
