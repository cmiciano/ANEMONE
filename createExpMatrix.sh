#!/bin/bash

### Take in a pwm matrix as argument one, genome as argument two
mpwm=$1
genome=$2

allmotifs="known"$genome".txt"

#Scan whole genome for motif occurrences

#scanMotifGenomeWide.pl $mpwm $genome -p 12 
scanMotifGenomeWide.pl $mpwm $genome -p 12 > $allmotifs 2> "known"$genome".log.txt"
#Separate the motif occurrences into separate files, annotate,
#and create a column of counts for each TF
#then create an expression matrix that counts the occurrence of each motif in the promoter
## of each gene
./sepMotifs.sh $allmotifs $genome
