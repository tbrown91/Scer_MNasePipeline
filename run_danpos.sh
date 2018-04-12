#!/bin/bash

python ./danpos-2.2.2/danpos.py dpos ./bam_files/BY_GLU/ -o ./bam_files/DANPOS_OUTPUT/BY_GLU/ -jd 146 -a 1
python ./danpos-2.2.2/danpos.py dpos ./bam_files/BY_15mGAL/ -o ./bam_files/DANPOS_OUTPUT/BY_15mGAL/ -jd 146 -a 1
python ./danpos-2.2.2/danpos.py dpos ./bam_files/BY_60mGAL/ -o ./bam_files/DANPOS_OUTPUT/BY_60mGAL/ -jd 146 -a 1


