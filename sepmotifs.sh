#!/bin/bash
file="$1"

## Script should take in the particularly genome, take the genomename and save it for later annotation

### Scanmotifs 
## scanMotifGenomeWide.pl known.motifs mm10 -p 12 > knownmm10.txt 2> knownmm10.log.txt

#awk -F/ '{print > $1".m"}' $file
./sepRAR.R 2> sepRAR.log.txt
#RAR\:RXR\(NR\)\,DR5\|ES-RAR-ChIP-Seq\(GSE56893\)-12
rm RAR\:RXR\(NR\),DR5\|ES-RAR-ChIP-Seq\(GSE56893\).m
echo "Finished separating RAR:RXR"
#awk -F/ '{print > $1"|"$2".m" }' $file
for tf in *.m; do
#for tf in *.txt; do
#for tf in anno*.txt; do
	#Adding anno prefix to TFname (anno(TF).m)
	#echo "Annotating" $tf
        antfm="anno"$tf

        #Output file name is anno(TF).txt
        antfname=${antfm%.m} #removing .m from endingi
	echo "Annotating" $tf
	fullanfile=$antfname".txt"
	echo "fullanfile" $fullanfile
        annotatePeaks.pl $tf mm10 > $fullanfile 2> $antfname".log.txt"


	## if the file is RAR:RXR, after annotation, separate the files        
	echo "Writing TF-countdf table for " $fullanfile
	./countMotifOccurCHIP.R $fullanfile 2> $antfname".count.log.txt"
	
#	echo "Parsing.." 
#        #./parsemotifs.R "annotzscan22nolabel.txt" > "parsestdout.txt" 2> "parsezscan.log.txt"
#        #./parsemotifs.R $suffix".txt" > $suffix".out.txt" 2> $suffix".log.txt"
done

#### Write the annotated file to a TF-countdf file
#for filename in anno*.txt; do
#        echo "Annotating" $filename
#        ./countMotifOccur.R $filename 2> $filename".log.txt"
#done

### Merge all of the TF-countdf files
echo "Merging files to make allTFcounts.txt"
./mergeFiles.R 2> "mergeFiles.log.txt"


