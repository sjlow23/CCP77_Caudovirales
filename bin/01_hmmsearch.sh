#!/bin/bash

# Script name: 01_hmmsearch.sh
# Description: Run hmmsearch and sselect top hits
# Note: Headers for protein fasta file must be proteinid followed by genomeid, separated by two dashes


usage() {
        echo -e "\nUsage: $0 [-i input_protein_fasta] [-j job_name]\n"
        echo 	"	-i      protein fasta file (required)"
        echo 	"	-j      job name (required)"
        echo -e "	-h      display help\n"
        exit 0
}

if [ ! -f 77markers.hmm ]
then
	echo -e "\nError: Can't find marker HMM file! Please check that file is in directory.\n"
	exit 0
fi

if [ $# -ne 4 ] 
then
        usage
        exit 0
else
        while getopts ":i:j:h" opt
        do
                case $opt in
                        i) IN=$OPTARG;;
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

mkdir 01_output
OUTDIR="01_output"

## Step 1- Run HMMsearch

echo -e "\nStep 1: Running HMMsearch..."
hmmsearch -o $OUTDIR/$IN.tsv --tblout 01_output/$IN.tbl --domtbl 01_output/$IN.domtbl --noali --notextw -E 1e-03 77markers.hmm $IN
echo -e "Step 1: HMMsearch complete!\n"


## Step 2- For identical genes hitting multiple markers, get top hit for each gene

echo -e "\nStep 2: Sorting best hits.."

cd $OUTDIR
grep -v "^#" $IN.domtbl | \
	awk '{ print $1, $4, $12, $14, ($17-$16+1)/$6*100, ($19-$18+1)/$3*100 }' OFS="\t" | \
	awk '{ if ($5 >= 50 && $6 >= 50) print $0 }' | \
	sort -k1,1 -k4rn | \
	awk '$1 != prev { print; prev = $1 }' > $IN.domtbl.filtered
sed 's/--/\t/1' $IN.domtbl.filtered | \
	awk -F "\t" '{ print $2"--"$3, $1, $4, $5, $6, $7 }' OFS="\t" | \
	sort -k1,1 -k4rn | \
	awk '$1 != prev { print; prev = $1 }' > $JOB.tab
sed 's/--/\t/1' $JOB.tab | \
	awk -F "\t" '{ print $3"--"$1, $2, $4, $5 }' OFS="\t" > $JOB.tab.final.txt
rm *.tsv *.tbl *.tab *.filtered
cd ..

echo -e "Step 2: Results table written to: 01_output/$job.tab.final.txt.\n"
