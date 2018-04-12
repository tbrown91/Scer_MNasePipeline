#!/bin/bash

bedtools genomecov -bga -ibam BY_D/BY_D_rep1_sorted.bam >& ../bedDraphs/BY_D_rep1.bed &
bedtools genomecov -bga -ibam BY_D/BY_D_rep2_sorted.bam >& ../bedDraphs/BY_D_rep2.bed &
bedtools genomecov -bga -ibam BY_15mG/BY_15mG_rep1_sorted.bam >& ../bed15mGraphs/BY_15mG_rep1.bed &
bedtools genomecov -bga -ibam BY_15mG/BY_15mG_rep2_sorted.bam >& ../bed15mGraphs/BY_15mG_rep2.bed &
bedtools genomecov -bga -ibam BY_60mG/BY_60mG_rep1_sorted.bam >& ../bed60mGraphs/BY_60mG_rep1.bed &
bedtools genomecov -bga -ibam BY_60mG/BY_60mG_rep2_sorted.bam >& ../bed60mGraphs/BY_60mG_rep2.bed &
