#!/bin/bash 

# Script name: 02_preprocess.sh
# Description: Extract sequences from HMMsearch output table for each marker and place in respective directories

#module load cdb/fasta emboss 

usage() {
        echo -e "\nUsage: $0 [-i HMM results tab file] [-p protein_fasta_file] [-m marker_list] [-j job_name]\n"
        echo 	"	-i      Results table from HMMsearch in step 1 (required)"
        echo 	"	-p      protein fasta file (required)"
        echo 	"	-m      list of selected markers (required)"
	echo 	"	-j	job name"
        echo -e "	-h      display help\n"
        exit 0
}


if [ $# -ne 8 ]
then
        usage
        exit 0
else
        while getopts ":i:p:m:j:h" opt
        do
                case $opt in
                        i) IN=$OPTARG;;
                        p) PROT=$OPTARG;;
                        m) MARKERS=$OPTARG;;
			j) JOB=$OPTARG;;
                        h) usage;;
                        \?) echo -e "\nInvalid option: -$OPTARG"
                        usage;;
                        :) echo -e "\nOption -$OPTARG requires an argument!\n"
                           exit 0;;
                esac
        done
        shift $(( OPTIND - 1 ))
fi


mkdir 02_output
OUTDIR="02_output"


index_and_extract() {
        cdbfasta $PROT
        cat $OUTDIR/$JOB.seq | cdbyank $PROT.cidx > $OUTDIR/"$PROT"_hits.faa
	rm $PROT.cidx
}

rename_fasta_headers() {
        scripts/translate_fasta_headers.pl --in=$OUTDIR/"$PROT"_hits.faa --tabfile=$OUTDIR/tabfile --out=$OUTDIR/"$JOB"_hits_final.faa
}


# Create tab file for renaming fasta headers

echo -e "\nStep 1- Renaming headers and extracting sequences.."

awk -F "\t" '{ print $1, $2 }' OFS="\t" $IN | \
	sed 's/--/\t/1' | \
	awk -F "\t" '{ print $1"--"$2, $3"--"$2 }' OFS="\t" > $OUTDIR/tabfile
awk '{ print $1 }' $PROT > $OUTDIR/$PROT.clean
mv $OUTDIR/$PROT.clean $OUTDIR/$PROT
awk '{ print $1 }' $IN > $OUTDIR/$JOB.seq
index_and_extract $PROT


# Rename headers of extracted sequences and remove unwanted files
rename_fasta_headers
echo -e "\nStep 1 - Done!\n"


# Create a new directory for each marker
cat $MARKERS | xargs -L1 mkdir
cdbfasta $OUTDIR/"$JOB"_hits_final.faa

# Loop each marker directory
echo -e "\nStep 2- Creating fasta files for each marker.."

while read vog
do
        grep $vog $OUTDIR/"$JOB"_hits_final.faa | sed 's/>//g' > $OUTDIR/$vog.seq
        cat $OUTDIR/$vog.seq | cdbyank $OUTDIR/"$JOB"_hits_final.faa.cidx > $vog.fasta
        mv $vog.fasta $vog
        
        cd $vog
        seqretsplit -sequence $vog.fasta -outseq $vog.out
        grep ">" $vog.fasta | \
		sed 's/>//g' | \
		awk '{ print $1".fasta" }' | sort > post
        ls *--*.fasta > pre
        paste pre post | xargs -L1 mv
        echo -e ""$vog": File renaming complete!\n"
        rm pre post
        cd ..

        rm $OUTDIR/$vog.seq
done < $MARKERS

echo -e "\nStep 2- Done!\n"
