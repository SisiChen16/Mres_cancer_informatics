#!/bin/bash
# Script that convert BAM files into fastq files and do the transcript isoform quantification using Salmon


## Usage Function
usage()
{
cat << EOF
usage: $0 options
This script is used to convert BAM files into a paired FastQ file for re-alignment to the transcriptome and read counting done by Salmon.

REQUIRED ARGUMENTS:
        -d 	Directory of where input BAM files are stored.

EXAMPLE USAGE:
   >> ~/02_Salmon_align.sh -d /path/to/data/dir/ 

EOF
}

## Parsing input arguments

DIR=
while getopts “hf:d:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         d)
             DIR=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done


## Handle empty arguments
if [[ -z $DIR ]]; then usage; exit 1; fi

cd $DIR 

for d in */ ; 
do
 echo "The CancerType processing is $d"
 INDIR=$DIR/$d
 cd $INDIR 
 for file in *.bam 
 	do 
 	OUTDIR=${file%%.*} 
 	mkdir -p $OUTDIR
 	echo "Start sorting BAM file $OUTDIR and converting them into fastq files"
	samtools sort -n $file -o $OUTDIR/sorted.bam
	samtools fastq -1 $OUTDIR/pair1.fq -2 $OUTDIR/pair2.fq -0 /dev/null -s /dev/null -n $OUTDIR/sorted.bam

	echo "Starting Salmon alignment on $OUTDIR on "`date`
	salmon quant -i /Users/chensisi/Documents/salmon_tutorial/Homo.sapiens_index -l A -1 $OUTDIR/pair1.fq -2 $OUTDIR/pair2.fq -p 8 --validateMappings -o $OUTDIR/salmon_quants
 	
 	echo "Deleting temporary files ..." 
 	rm $OUTDIR/pair1.fq
 	rm $OUTDIR/pair2.fq
 	rm $OUTDIR/sorted.bam

 		echo "Finished sequence alignment of $OUTDIR on "`date`
 done

done

echo "Completed Isoform quantification on "`date`









