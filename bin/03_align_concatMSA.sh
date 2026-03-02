#!/bin/bash 

# Script name: run_hmmalign.sh
# Description: Aligns protein sequences to their VOG profile hmms and converts to string format
# Description: Concatenate all single marker alignments, trim and generate final super alignment


#module load hmmer cdb/fasta emboss parallel biosquid trimal

usage() {
        echo -e "\nUsage: $0 [-m marker list] [-t marker count threshold] [-j job name]\n"
        echo 	"	-m      list of markers (required)"
	echo	"	-t	minimum number of markers present in genome to pass filtering (required)"
	echo 	"	-j	job name"
        echo -e "	-h      display help\n"
        exit 0
}


if [ $# -ne 6 ]
then
        usage
        exit 0
else
        while getopts ":m:t:j:h" opt
        do
                case $opt in
                        m) MARKERS=$OPTARG;;
			t) COUNT=$OPTARG;;
			j) JOB=$OPTARG;;
                        h) usage;;
                        :) echo -e "\nOption -$OPTARG requires an argument!\n"
                           exit 0;;
                        \?) echo -e "\nInvalid option: -$OPTARG"
                            usage;;
                esac
        done
        shift $(( OPTIND - 1 ))
fi

mkdir 03_output
OUTDIR="03_output"


align_to_hmm() {
        echo ""$vog": hmmalign in progress..."
        ls *--*.fasta | awk '{ print "hmmalign --outformat Pfam '"$vog"'.hmm "$0" > "$0".hmmaln" }' > align.sh
        cat align.sh | parallel -j 20 
        rm *--*.fasta
        rm align.sh
        echo -e ""$vog": hmmalign complete!\n"
}

hmmaln_to_string() {
        ls *.hmmaln | awk '{ print "../scripts/hmmaln_to_string.py "$0, $0" > "$0".string" }' | \
		sed 's/.fasta.hmmaln//1' | sed 's/fasta.hmmaln.//2' > aln2string.sh
        echo -e ""$vog": Converting to string..\n"
        cat aln2string.sh | parallel -j 8

        for i in *.string; do awk -F, 'NR==0{print "file_name",$0;next}{print FILENAME , $0}' OFS="\t" $i | sed 's/\.string//g' > $i.tmp; done
        cat VOG*.tmp | awk '{ print ">"$0"" }' | \
		sed 's/\t/\n/g' > "$vog".hmmaln.faa
        rm *.hmmaln aln2string.sh VOG*.string *.tmp
}

trim_aln_t50() {
        echo -e ""$vog": Trimming alignment...\n"
        trimal -in $vog.hmmaln.faa -out $vog.hmmaln.faa.trim -gt 0.5
}

missing_markers() {
	echo -e ""$vog": Generating gapped sequences...\n"
        sed -i "s/"$vog"--//g" *.trim
        echo "-" > tmp
        length="`seqstat *.trim | grep Small | awk -F ":" '{ print $2 }' | sed 's/ //g'`"
        perl -ne 'print "$_" x'$length'' tmp | \
		tr -d "\n" | awk '{ print $0"\t" }' > gaps
        grep ">" *.trim | awk '{ print $1 }' | \
		sed 's/>//g' | sort > $vog.headers.sorted
        comm -23 ../$OUTDIR/genomes_min$COUNT $vog.headers.sorted > tofill
        rm tmp $vog.headers.sorted
}

process() {
        while read line
        do
                echo -e ">"$line"\n$(cat gaps)" >> $vog.hmmaln.faa.trim
        done < tofill
        rm tofill gaps
}

get_aln() {
        cdbfasta $vog.hmmaln.faa.trim
        cat ../$OUTDIR/genomes_min$COUNT | \
		cdbyank $vog.hmmaln.faa.trim.cidx > $vog.hmmaln.faa.trim.sorted
        ../scripts/fasta_MLtoSL.py $vog.hmmaln.faa.trim.sorted
	
	#awk '{ print $1 }' $vog.hmmaln.faa.trim.sorted.out > $vog.hmmaln.faa.trim.sorted.clean
	#mv $vog.hmmaln.faa.trim.sorted.clean $vog.hmmaln.faa.trim
        grep -v ">" $vog.hmmaln.faa.trim > col2
        rm *.cidx *.sorted
}


##Untar 77markers.tar.gz and move to respective directories
tar -xzvf 77markers.tar.gz
ls -d VOG*.hmm | awk '{ print $1, $1 }' | sed 's/.hmm//2' | xargs -L1 mv


##Align and trim here
while read vog
do
	cd $vog
	align_to_hmm
	hmmaln_to_string
	trim_aln_t50
	cd ..
done < $MARKERS


grep ">" VOG*/*.trim | \
	awk -F ">" '{ print $2 }' | awk '{ print $1 }' | \
	awk -F "--" '{ print $2 }' | \
	sort | uniq -c | awk '{ if ($1 >= '$COUNT') print $2 }' > $OUTDIR/genomes_min$COUNT
echo -e "\nAll alignments complete!\n"

##Generate gapped sequences and concatenate individual marker MSAs here
while read vog
do
        cd $vog
        missing_markers
        process
        get_aln
        cd ..
done < $MARKERS


paste VOG*/col2 | sed 's/\t//g' > $OUTDIR/$JOB.col2.aln
rm VOG*/col2

##Concatenate and trim final alignment
cd $OUTDIR
	cat genomes_min$COUNT | \
		awk '{ print ">"$0"" }' > $JOB.col1.aln
	paste $JOB.col1.aln $JOB.col2.aln | \
		sed 's/\t/\n/g' > CCP77.aln
	msa_length="`seqstat CCP77.aln | grep Small | awk -F ":" '{ print $2 }' | sed 's/ //g'`"

	echo -e "\nGetting final alignment."
	degapseq -sequence CCP77.aln -outseq CCP77.aln.degap
	seqstat -a CCP77.aln.degap | \
		grep "*" | awk '{ if ($3 >= 0.05*'$msa_length') print $2 }' > genomes_5p
	cdbfasta CCP77.aln && cat genomes_5p | cdbyank CCP77.aln.cidx > CCP77.trimmed.aln
	rm $JOB.*.aln CCP77.aln.cidx $OUTDIR/CCP77.aln.degap
cd ..


##Gzip all marker folders
tar --remove-files -czvf VOGalignments.tar.gz VOG*/
echo -e "\nFinal alignment ready."
